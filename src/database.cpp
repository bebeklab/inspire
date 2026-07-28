#include "database.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
class StatementGuard {
public:
    explicit StatementGuard(sqlite3_stmt* stmt) : stmt_(stmt) {}
    ~StatementGuard() { if (stmt_) sqlite3_finalize(stmt_); }
    sqlite3_stmt* get() const { return stmt_; }
private:
    sqlite3_stmt* stmt_;
};

void bindText(sqlite3_stmt* stmt, int index, const std::string& value) {
    const int rc = sqlite3_bind_text(
        stmt, index, value.c_str(), -1, SQLITE_TRANSIENT);
    if (rc != SQLITE_OK) {
        throw std::runtime_error("sqlite3_bind_text failed");
    }
}
}

InspireDatabase::InspireDatabase(const std::string& fileName) : db_(NULL) {
    const int flags = SQLITE_OPEN_READONLY | SQLITE_OPEN_FULLMUTEX;
    const int rc = sqlite3_open_v2(fileName.c_str(), &db_, flags, NULL);
    if (rc != SQLITE_OK) {
        const std::string message = db_ ? sqlite3_errmsg(db_) : "unknown error";
        if (db_) sqlite3_close(db_);
        db_ = NULL;
        throw std::runtime_error("Cannot open SQLite database: " + message);
    }
    sqlite3_busy_timeout(db_, 30000);
    sqlite3_extended_result_codes(db_, 1);
}

InspireDatabase::~InspireDatabase() {
    if (db_) sqlite3_close(db_);
}

[[noreturn]] void InspireDatabase::throwError(
    const std::string& context) const {
    throw std::runtime_error(context + ": " + sqlite3_errmsg(db_));
}

sqlite3_stmt* InspireDatabase::prepare(const char* sql) const {
    sqlite3_stmt* stmt = NULL;
    const int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, NULL);
    if (rc != SQLITE_OK) throwError("Preparing SQLite statement");
    return stmt;
}

std::string InspireDatabase::columnText(sqlite3_stmt* stmt, int column) {
    const unsigned char* text = sqlite3_column_text(stmt, column);
    return text ? reinterpret_cast<const char*>(text) : std::string();
}

std::vector<EdgeRow> InspireDatabase::loadEdges(
    const std::string& petal) const {
    StatementGuard guard(prepare(
        "SELECT gene1, gene2, interaction_type, max_score "
        "FROM edge WHERE petal=?1 ORDER BY gene1, gene2"));
    bindText(guard.get(), 1, petal);

    std::vector<EdgeRow> rows;
    int rc = SQLITE_ROW;
    while ((rc = sqlite3_step(guard.get())) == SQLITE_ROW) {
        EdgeRow row;
        row.gene1 = columnText(guard.get(), 0);
        row.gene2 = columnText(guard.get(), 1);
        row.interactionType = columnText(guard.get(), 2);
        row.ppiScore = sqlite3_column_type(guard.get(), 3) == SQLITE_NULL
            ? 0.0 : sqlite3_column_double(guard.get(), 3);
        rows.push_back(row);
    }
    if (rc != SQLITE_DONE) throwError("Loading petal edges");
    return rows;
}

std::vector<std::string> InspireDatabase::loadPetalGenes(
    const std::string& petal) const {
    StatementGuard guard(prepare(
        "SELECT gene_symbol FROM petal_gene "
        "WHERE petal=?1 ORDER BY gene_symbol"));
    bindText(guard.get(), 1, petal);

    std::vector<std::string> genes;
    int rc = SQLITE_ROW;
    while ((rc = sqlite3_step(guard.get())) == SQLITE_ROW) {
        genes.push_back(columnText(guard.get(), 0));
    }
    if (rc != SQLITE_DONE) throwError("Loading petal genes");
    return genes;
}

std::string InspireDatabase::preferredUniProt(
    const std::string& geneSymbol) const {
    StatementGuard guard(prepare(
        "SELECT uniprot_accession FROM gene_uniprot "
        "WHERE gene_symbol=?1 AND is_preferred=1"));
    bindText(guard.get(), 1, geneSymbol);

    const int rc = sqlite3_step(guard.get());
    if (rc == SQLITE_DONE) {
        throw std::runtime_error(
            "No preferred UniProt accession for gene " + geneSymbol);
    }
    if (rc != SQLITE_ROW) throwError("Loading preferred UniProt accession");
    const std::string value = columnText(guard.get(), 0);
    if (sqlite3_step(guard.get()) == SQLITE_ROW) {
        throw std::runtime_error(
            "Multiple preferred UniProt accessions for gene " + geneSymbol);
    }
    return value;
}

std::vector<std::string> InspireDatabase::loadPfam(
    const std::string& uniprot) const {
    StatementGuard guard(prepare(
        "SELECT pfam_accession FROM protein_pfam "
        "WHERE uniprot_accession=?1 ORDER BY pfam_accession"));
    bindText(guard.get(), 1, uniprot);

    std::vector<std::string> pfams;
    int rc = SQLITE_ROW;
    while ((rc = sqlite3_step(guard.get())) == SQLITE_ROW) {
        pfams.push_back(columnText(guard.get(), 0));
    }
    if (rc != SQLITE_DONE) throwError("Loading Pfam annotations");
    return pfams;
}

double InspireDatabase::loadGOScore(
    const std::string& uniprotA,
    const std::string& uniprotB,
    const std::string& algorithm) const {
    const std::string first = std::min(uniprotA, uniprotB);
    const std::string second = std::max(uniprotA, uniprotB);

    StatementGuard guard(prepare(
        "SELECT score FROM go_pair_score "
        "WHERE protein1=?1 AND protein2=?2 AND algorithm=?3 "
        "ORDER BY rowid DESC LIMIT 1"));
    bindText(guard.get(), 1, first);
    bindText(guard.get(), 2, second);
    bindText(guard.get(), 3, algorithm);

    const int rc = sqlite3_step(guard.get());
    if (rc == SQLITE_DONE) {
        throw std::runtime_error(
            "Missing GO pair score for " + first + ", " + second +
            " using " + algorithm);
    }
    if (rc != SQLITE_ROW) throwError("Loading GO pair score");
    return sqlite3_column_double(guard.get(), 0);
}

double InspireDatabase::maximumGOScore(const std::string& algorithm) const {
    StatementGuard guard(prepare(
        "SELECT MAX(score) FROM go_pair_score WHERE algorithm=?1"));
    bindText(guard.get(), 1, algorithm);
    const int rc = sqlite3_step(guard.get());
    if (rc != SQLITE_ROW) throwError("Loading GO normalization value");
    const double value = sqlite3_column_double(guard.get(), 0);
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::runtime_error("Invalid maximum GO score for " + algorithm);
    }
    return value;
}

int InspireDatabase::preferredProteinCount() const {
    StatementGuard guard(prepare(
        "SELECT COUNT(DISTINCT uniprot_accession) "
        "FROM gene_uniprot WHERE is_preferred=1"));
    if (sqlite3_step(guard.get()) != SQLITE_ROW) {
        throwError("Counting preferred proteins");
    }
    return sqlite3_column_int(guard.get(), 0);
}

bool InspireDatabase::hasCompleteGOScoreCache(
    const std::string& algorithm) const {
    const long long n = preferredProteinCount();
    const long long expected = n * (n + 1) / 2;

    StatementGuard guard(prepare(
        "SELECT COUNT(*) FROM go_pair_score WHERE algorithm=?1"));
    bindText(guard.get(), 1, algorithm);
    if (sqlite3_step(guard.get()) != SQLITE_ROW) {
        throwError("Counting GO scores");
    }
    const long long actual = sqlite3_column_int64(guard.get(), 0);
    return actual == expected;
}

std::string InspireDatabase::metadata(const std::string& key) const {
    StatementGuard guard(prepare(
        "SELECT value FROM metadata WHERE key=?1"));
    bindText(guard.get(), 1, key);
    const int rc = sqlite3_step(guard.get());
    if (rc == SQLITE_DONE) return std::string();
    if (rc != SQLITE_ROW) throwError("Loading metadata");
    return columnText(guard.get(), 0);
}
