#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/IMeteorite.h"

// Collection with information about popular meteorites
// Implemented as the singleton class
class KnownMeteorites
{
  public:
    // IDs of several meteorites that have been fairly well documented and studied
    // To access other records, use string-based identifiers
    enum class ID
    {
      PRIBRAM,
      LOST_CITY,
      INNISFREE,
      PEEKSKILL,
      KOSICE,
      CHELYABINSK
    };

    KnownMeteorites(const KnownMeteorites &) = delete;
    KnownMeteorites &operator =(const KnownMeteorites &) = delete;

    // Simplified version that returns information about a predefined meteorite
    // Feel free to use it for testing and debugging purposes
    const std::shared_ptr<const IMeteorite> &Get(ID id) const;

    // Queries all records about meteorites with the specified name AND the required DOI
    // Use an empty string if you accept any
    // In some cases, different records can describe the same meteorite
    // This means that the used information is taken from different sources (tables, DOIs, etc)
    const std::vector<std::shared_ptr<const IMeteorite>> Get(std::string name = std::string(),
                                                         std::string doi = std::string()) const;

    // Returns names of all known meteorites (names are unique)
    const std::vector<std::string> &Names() const
    { return names_; }

    // Returns DOIs of all research papers from which we extracted information about known meteorites
    const std::vector<std::string> &DOIs() const
    { return dois_; }

    // Provides the reference to the singleton instance of this class
    static const KnownMeteorites &Ref();

  private:
    KnownMeteorites();

    std::unordered_map<ID, std::shared_ptr<const IMeteorite>> predefined_;
    std::vector<std::string> names_, dois_;
    std::vector<std::shared_ptr<const IMeteorite>> collection_;
};
