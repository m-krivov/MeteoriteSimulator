#pragma once
#include "Meteorites.Core/Defs.h"

#include "Meteorites.Core/IMeteorite.h"

class KnownMeteorites
{
  public:
    enum class ID : uint32_t
    {
      PRIBRAM   = 0,
      LOST_CITY = 1,
      INNISFREE = 2,
      KOSICE    = 3
    };

    KnownMeteorites() = delete;
    KnownMeteorites(const KnownMeteorites &) = delete;
    KnownMeteorites &operator =(const KnownMeteorites &) = delete;

    static const IMeteorite *Get(ID id);

    static const std::vector<const IMeteorite *> &All();

  private:
    static struct Collection
    {
      std::unordered_map<ID, std::unique_ptr<IMeteorite> > meteorites;
      std::vector<const IMeteorite *> enumeration;
    } collection_;

    static void LazyInit();
};
