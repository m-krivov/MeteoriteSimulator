#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/IMeteorite.h"

class KnownMeteorites
{
  public:
    enum MeteoriteID
    {
      INNISFREE = 0
    };

    KnownMeteorites() = delete;
    KnownMeteorites(const KnownMeteorites &) = delete;
    void operator =(const KnownMeteorites &) = delete;

    static const IMeteorite *Get(MeteoriteID id);

    static const std::vector<const IMeteorite *> &All();

  private:
    static struct Collection
    {
      std::unordered_map<MeteoriteID, std::unique_ptr<IMeteorite> > meteorites;
      std::vector<const IMeteorite *> enumeration;
    } collection_;

    static void LazyInit();
};
