#ifndef INSPIRE_DATABASE_H
#define INSPIRE_DATABASE_H

#include <sqlite3.h>

#include <string>
#include <utility>
#include <vector>

struct EdgeRow {
    std::string gene1;
    std::string gene2;
    std::string interactionType;
    double ppiScore;
};

class InspireDatabase {
public:
    explicit InspireDatabase(const std::string& fileName);
    ~InspireDatabase();

    InspireDatabase(const InspireDatabase&) = delete;
    InspireDatabase& operator=(const InspireDatabase&) = delete;

    std::vector<EdgeRow> loadEdges(const std::string& petal) const;
    std::vector<std::string> loadPetalGenes(const std::string& petal) const;
    std::string preferredUniProt(const std::string& geneSymbol) const;
    std::vector<std::string> loadPfam(const std::string& uniprot) const;

    double loadGOScore(
        const std::string& uniprotA,
        const std::string& uniprotB,
        const std::string& algorithm) const;

    double maximumGOScore(const std::string& algorithm) const;
    bool hasCompleteGOScoreCache(const std::string& algorithm) const;
    int preferredProteinCount() const;
    std::string metadata(const std::string& key) const;

private:
    sqlite3* db_;

    sqlite3_stmt* prepare(const char* sql) const;
    static std::string columnText(sqlite3_stmt* stmt, int column);
    [[noreturn]] void throwError(const std::string& context) const;
};

#endif
