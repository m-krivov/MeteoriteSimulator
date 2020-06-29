#include "KnownMeteorites.h"

namespace
{

// Basic class that implements interface 'IMeteorites' and overrides its optional fields by default values
class BasicMeteorite : public IMeteorite
{
  public:
    virtual void Trajectory(size_t &records, const real *&time,
                            const real *&v, const real *&h) const override
    {
      if (time_.empty())
      {throw std::runtime_error("Internal error, class 'BasicMeteorite' was not configured"); }

      records = time_.size();
      time = &time_[0];
      v = &v_[0];
      h = &h_[0];
    }

  protected:
    BasicMeteorite() = default;

    void SetTrajectory(std::vector<real> &time, std::vector<real> &v, std::vector<real> &h)
    {
      if (time.size() != v.size() || time.size() != h.size())
      { throw std::runtime_error("Tables for time, velocity and height must have the same size"); }
      if (time.empty())
      { throw std::runtime_error("Tables for time, velocity and height must have at least one record"); }

      time_ = std::move(time);
      v_ = std::move(v);
      h_ = std::move(h);
    }

  private:
    std::vector<real> time_, v_, h_;
};

}

//-----------------
//--- Innisfree ---
//-----------------

namespace
{

class InnisfreeMeteorite : public BasicMeteorite
{
  public:
    virtual std::string Name() const override
    { return "Innisfree"; }

    virtual std::string Date() const override
    { return "5th of February, 1977"; }

    virtual std::string FallLocation() const override
    { return "Innisfree, Alberta, Canada"; }

    InnisfreeMeteorite()
    {
      // 1981_Halliday
      std::vector<real> time = { (real)0.0, (real)0.2, (real)0.4, (real)0.6, (real)0.8, (real)1.0,
                                 (real)1.2, (real)1.4, (real)1.6, (real)1.8, (real)2.0, (real)2.2,
                                 (real)2.4, (real)2.6, (real)2.8, (real)3.0, (real)3.2 };
      std::vector<real> h    = { (real)58800, (real)56100, (real)53500, (real)50800, (real)48200, (real)45500,
                                 (real)42800, (real)40200, (real)37500, (real)35000, (real)32500, (real)30200,
                                 (real)27900, (real)25900, (real)24200, (real)22600, (real)21500 };
      std::vector<real> v    = { (real)14510, (real)14490, (real)14470, (real)14440, (real)14340, (real)14230,
                                 (real)14050, (real)13790, (real)13420, (real)12960, (real)12350, (real)11540,
                                 (real)10430, (real)8890,  (real)7240,  (real)5540,  (real)4700 };
      SetTrajectory(time, v, h);
    }
};

}

//-----------------------
//--- KnownMeteorites ---
//-----------------------

KnownMeteorites::Collection KnownMeteorites::collection_;

void KnownMeteorites::LazyInit()
{
  if (collection_.enumeration.empty())
  {
    std::unique_ptr<IMeteorite> innisfree(new InnisfreeMeteorite());
    collection_.enumeration.emplace_back(innisfree.get());
    collection_.meteorites.emplace(std::make_pair(INNISFREE, std::move(innisfree)));
  }
}

const IMeteorite *KnownMeteorites::Get(MeteoriteID id)
{
  LazyInit();
  auto iter = collection_.meteorites.find(id);
  if (iter == collection_.meteorites.end())
  {
    std::stringstream ss;
    ss << "Internal error, meteorite with ID=" << id << " is not found";
    throw std::runtime_error(ss.str());
  }

  return iter->second.get();
}

const std::vector<const IMeteorite *> &KnownMeteorites::All()
{
  LazyInit();
  return collection_.enumeration;
}
