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
      { throw std::runtime_error("Internal error, class 'BasicMeteorite' was not configured"); }

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

//---------------
//--- Pribram ---
//---------------

namespace
{

class PribramMeteorite : public BasicMeteorite
{
  public:
    virtual std::string Name() const override
    { return "Pribram"; }

    virtual std::string Date() const override
    {return "7th of April, 1959"; }

    virtual std::string FallLocation() const override
    { return " Pribram, Czechoslovakia (Czech Republic)"; }

    PribramMeteorite()
    {
      // 2008_Gritsevich.pdf
      std::vector<real> time = { (real)0,      (real)0.856,  (real)1.732,
                                 (real)2.494,  (real)2.692,  (real)3.0676 };

      std::vector<real> h    = { (real)88.594, (real)76.318, (real)63.837,
                                 (real)52.97,  (real)50.164, (real)44.858 };
      for (auto &r : h) { r *= 1000; }

      std::vector<real> v    = { (real)20.887, (real)20.86,  (real)20.838,
                                 (real)20.773, (real)20.717, (real)20.459 };
      for (auto &r : v) { r *= 1000; }

      SetTrajectory(time, v, h);
    }
};

}

namespace
{

class LostCityMeteorite : public BasicMeteorite
{
  public:
    virtual std::string Name() const override
    { return "Lost City"; }

    virtual std::string Date() const override
    { return "3rd of January, 1970"; }

    virtual std::string FallLocation() const override
    { return "Lost City, Oklahoma, USA"; }

    LostCityMeteorite()
    {
      // 2008_Gritsevich.pdf
      std::vector<real> time = { (real)0.05, (real)1.05, (real)2.05, (real)3.05, (real)4.05,
                                 (real)5.05, (real)6.05, (real)7.05, (real)8.00, (real)8.95 };
      {
        auto t0 = time[0];
        for (auto &r : time) { r -= t0; }
      }

      std::vector<real> h    = { (real)85.9, (real)77.1, (real)68.5, (real)59.9, (real)51.3,
                                 (real)42.8, (real)34.6, (real)27.5, (real)22.6, (real)19.9 };
      for (auto &r : h) { r *= 1000; }

      std::vector<real> v    = { (real)14.2, (real)14.2, (real)14.2, (real)14.1, (real)14.0,
                                 (real)13.8, (real)12.9, (real)10.3, (real)6.1,  (real)3.4 };
      for (auto &r : v) { r *= 1000; }

      SetTrajectory(time, v, h);
    }
};

} // unnamed namespace

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

} // unnamed namespace

//--------------
//--- Kosice ---
//--------------

namespace
{

class KosiceMeteorite : public BasicMeteorite
{
  public:
    virtual std::string Name() const override
    { return "Kosice"; }

    virtual std::string Date() const override
    { return "28th of February, 2010"; }

    virtual std::string FallLocation() const override
    { return "Kosice Region, Slovakia"; }

    KosiceMeteorite()
    {
      // 2017_Gritsevich, but sorted!
      std::vector<real> time = { (real)2.16,    (real)2.24,    (real)2.32,    (real)2.4,     (real)2.48,
                                 (real)2.56,    (real)2.64,    (real)2.72,    (real)2.8,     (real)2.88,
                                 (real)2.96,    (real)3.04,    (real)3.12,    (real)3.16169, (real)3.2,
                                 (real)3.28,    (real)3.36,    (real)3.36169, (real)3.44,    (real)3.52,
                                 (real)3.56169, (real)3.6,     (real)3.68,    (real)3.76169, (real)3.96169,
                                 (real)4.16169, (real)4.36169, (real)4.4,     (real)4.76169, (real)4.96169,
                                 (real)5.04,    (real)5.16169, (real)5.12,    (real)5.2,     (real)5.28,
                                 (real)5.36,    (real)5.36169, (real)5.44,    (real)5.52,    (real)5.6,
                                 (real)5.56169, (real)5.68,    (real)5.71304, (real)5.76,    (real)5.84,
                                 (real)5.91304, (real)5.92,    (real)6.0,     (real)6.08,    (real)6.11304,
                                 (real)6.16,    (real)6.24,    (real)6.31304, (real)6.32,    (real)6.4,
                                 (real)6.48,    (real)6.51304, (real)6.56,    (real)6.64,    (real)6.71304 };
      {
        auto t0 = time[0];
        for (auto &r : time) { r -= t0; }
      }

      std::vector<real> h =    { (real)68.2764, (real)68.0288, (real)66.5686, (real)64.5812, (real)63.9505,
                                 (real)63.7422, (real)62.6369, (real)61.4596, (real)60.4398, (real)59.4206,
                                 (real)58.4696, (real)57.3074, (real)56.2359, (real)55.1798, (real)54.9518,
                                 (real)54.6753, (real)52.7791, (real)52.2932, (real)51.9038, (real)51.2779,
                                 (real)50.4611, (real)49.6837, (real)49.4245, (real)47.6025, (real)45.1692,
                                 (real)42.4232, (real)39.9047, (real)39.0502, (real)35.2954, (real)32.9749,
                                 (real)30.544,  (real)30.4843, (real)30.2927, (real)29.7595, (real)29.2252,
                                 (real)27.817,  (real)27.3005, (real)26.4961, (real)25.6557, (real)25.4808,
                                 (real)25.4054, (real)23.78,   (real)22.6837, (real)22.6732, (real)22.2856,
                                 (real)21.7107, (real)21.3684, (real)21.0391, (real)20.6797, (real)20.5059,
                                 (real)19.9564, (real)19.4099, (real)19.2934, (real)19.1465, (real)18.8756,
                                 (real)18.6043, (real)18.3224, (real)18.3128, (real)18.0459, (real)17.4214 };
      for (auto &r : h) { r *= 1000; }

      std::vector<real> v =    { (real)14.97917, (real)14.97882, (real)14.97841, (real)14.97792, (real)14.97734,
                                 (real)14.97667, (real)14.97588, (real)14.97495, (real)14.973,   (real)14.97258,
                                 (real)14.97108, (real)14.96931, (real)14.96724, (real)14.96602, (real)14.9648,
                                 (real)14.96194, (real)14.9585,  (real)14.95858, (real)14.95463, (real)14.95,
                                 (real)14.94727, (real)14.94456, (real)14.93816, (real)14.93048, (real)14.90538,
                                 (real)14.86787, (real)14.81179, (real)14.79823, (real)14.60266, (real)14.41535,
                                 (real)14.31889, (real)14.20335, (real)14.13536, (real)14.06766, (real)13.90829,
                                 (real)13.72112, (real)13.71683, (real)13.5013,  (real)13.24313, (real)13.09121,
                                 (real)12.93992, (real)12.58382, (real)12.41921, (real)12.1656,  (real)11.67442,
                                 (real)11.1515,  (real)11.09755, (real)10.42005, (real)9.127726, (real)8.39702,
                                 (real)7.511363, (real)6.337642, (real)5.549147, (real)5.485346, (real)4.866452,
                                 (real)4.450264, (real)4.393314, (real)4.326966, (real)4.244316, (real)4.192904 };
      for (auto &r : v) { r *= 1000; }

      SetTrajectory(time, v, h);
    }
};

} // unnamed namespace

//-----------------------
//--- KnownMeteorites ---
//-----------------------

KnownMeteorites::Collection KnownMeteorites::collection_;

void KnownMeteorites::LazyInit()
{
  if (collection_.enumeration.empty())
  {
    std::unique_ptr<IMeteorite> pribram(new PribramMeteorite());
    collection_.enumeration.emplace_back(pribram.get());
    collection_.meteorites.emplace(std::make_pair(ID::PRIBRAM, std::move(pribram)));

    std::unique_ptr<IMeteorite> lostcity(new LostCityMeteorite());
    collection_.enumeration.emplace_back(lostcity.get());
    collection_.meteorites.emplace(std::make_pair(ID::LOST_CITY, std::move(lostcity)));

    std::unique_ptr<IMeteorite> innisfree(new InnisfreeMeteorite());
    collection_.enumeration.emplace_back(innisfree.get());
    collection_.meteorites.emplace(std::make_pair(ID::INNISFREE, std::move(innisfree)));

    std::unique_ptr<IMeteorite> kosice(new KosiceMeteorite());
    collection_.enumeration.emplace_back(kosice.get());
    collection_.meteorites.emplace(std::make_pair(ID::KOSICE, std::move(kosice)));
  }
}

const IMeteorite *KnownMeteorites::Get(ID id)
{
  LazyInit();
  auto iter = collection_.meteorites.find(id);
  if (iter == collection_.meteorites.end())
  {
    std::stringstream ss;
    ss << "Internal error, meteorite with ID=" << (uint32_t)id << " is not found";
    throw std::runtime_error(ss.str());
  }

  return iter->second.get();
}

const std::vector<const IMeteorite *> &KnownMeteorites::All()
{
  LazyInit();
  return collection_.enumeration;
}
