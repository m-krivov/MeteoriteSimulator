#include "KnownMeteorites.h"

#include "KnownMeteorites/1981_Halliday.h"
#include "KnownMeteorites/1995_Beech.h"
#include "KnownMeteorites/2008_Gritsevich.h"
#include "KnownMeteorites/2013_Borovicka.h"
#include "KnownMeteorites/2017_Gritsevich.h"

KnownMeteorites::KnownMeteorites()
{
  // Populate the collection with known meteorites
  Halliday1981::Populate(collection_);
  Beech1995::Populate(collection_);
  Gritsevich2008::Populate(collection_);
  Borovicka2013::Populate(collection_);
  Gritsevich2017::Populate(collection_);
  /* add here any new records */

  // Merge the names and DOIs of all meteorites
  std::unordered_set<std::string> names, dois;
  for (const auto &meteorite : collection_)
  {
    names.insert(meteorite->Name());
    dois.insert(meteorite->DOI());
  }
  names_ = std::move(std::vector<std::string>(names.begin(), names.end()));
  dois_ = std::move(std::vector<std::string>(dois.begin(), dois.end()));

  // Register the predefined meteorites by their identifiers
  {
    auto records = Get("Innisfree", "https://doi.org/10.1111/j.1945-5100.1981.tb00540.x");
    assert(records.size() >= 1);
    predefined_.emplace(ID::INNISFREE, records[0]);
  }
  {
    auto records = Get("Lost City", "https://doi.org/10.1134/S003809460805002X");
    assert(records.size() >= 1);
    predefined_.emplace(ID::LOST_CITY, records[0]);
  }
  {
    auto records = Get("Pribram", "https://doi.org/10.1134/S003809460805002X");
    assert(records.size() >= 1);
    predefined_.emplace(ID::PRIBRAM, records[0]);
  }
  {
    auto records = Get("Peekskill", "https://doi.org/10.1007/BF00671508");
    assert(records.size() >= 1);
    predefined_.emplace(ID::PEEKSKILL, records[0]);
  }
  {
    auto records = Get("Kosice", "https://doi.org/10.1007/978-3-319-46179-3_8");
    assert(records.size() >= 1);
    predefined_.emplace(ID::KOSICE, records[0]);
  }
  {
    auto records = Get("Chelyabinsk", "https://doi.org/10.1038/nature12671");
    assert(records.size() >= 1);
    predefined_.emplace(ID::CHELYABINSK, records[0]);
  }
}

const IMeteorite &KnownMeteorites::Get(ID id) const
{
  auto iter = predefined_.find(id);
  assert(iter != predefined_.end());
  assert(iter->second != nullptr);
  return *iter->second;
}

const std::vector<const IMeteorite *> KnownMeteorites::Get(std::string name, std::string doi) const
{
  std::vector<const IMeteorite *> result;
  for (const auto &meteorite : collection_)
  {
    bool name_is_ok = name.empty() || meteorite->Name() == name;
    bool doi_is_ok = doi.empty() || meteorite->DOI() == doi;
    if (name_is_ok && doi_is_ok)
    { result.emplace_back(meteorite.get()); }
  }
  return result;
}

const KnownMeteorites &KnownMeteorites::Ref()
{
  static std::shared_ptr<KnownMeteorites> instance;
  if (instance == nullptr)
  { instance.reset(new KnownMeteorites()); }
  return *instance;
}
