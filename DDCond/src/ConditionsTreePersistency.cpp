//==========================================================================
//  AIDA Detector description implementation 
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================

// Framework include files
#include <DD4hep/Printout.h>
#include <DDCond/ConditionsSlice.h>
#include <DDCond/ConditionsIOVPool.h>
#include <DDCond/ConditionsTreePersistency.h>

#include <TFile.h>
#include <TTimeStamp.h>

typedef dd4hep::cond::ConditionsTreePersistency __ConditionsTreePersistency;
/// ROOT Class implementation directive
ClassImp(__ConditionsTreePersistency)

using namespace dd4hep::cond;

// Local namespace for anonymous stuff
namespace  {

  /// Helper to select conditions
  /*
   *  \author  M.Frank
   *  \version 1.0
   *  \ingroup DD4HEP_CONDITIONS
   */
  struct Scanner : public dd4hep::Condition::Processor   {
    ConditionsTreePersistency::pool_type& pool;
    /// Constructor
    Scanner(ConditionsTreePersistency::pool_type& p) : pool(p) {}
    /// Conditions callback for object processing
    virtual int process(dd4hep::Condition c)  const override  {
      pool.emplace_back(c.ptr());
      return 1;
    }
  };
  struct DurationStamp  {
    TTimeStamp start;
    ConditionsTreePersistency* object = 0;
    DurationStamp(ConditionsTreePersistency* obj) : object(obj)  {
    }
    ~DurationStamp()  {
      TTimeStamp stop;
      object->duration = stop.AsDouble()-start.AsDouble();
    }
  };
}

/// Default constructor
ConditionsTreePersistency::ConditionsTreePersistency() : TNamed()   {
}

/// Initializing constructor
ConditionsTreePersistency::ConditionsTreePersistency(const std::string& name, const std::string& title)
  : TNamed(name.c_str(), title.c_str())
{
}

/// Default destructor
ConditionsTreePersistency::~ConditionsTreePersistency()    {
  clear();
}

/// Add conditions content to be saved. Note, that dependent conditions shall not be saved!
std::size_t ConditionsTreePersistency::add(const std::string& identifier,
                                           const IOV& iov,
                                           std::vector<Condition>& conditions)
{
  DurationStamp stamp(this);
  conditionPools.emplace_back(std::pair<iov_key_type, pool_type>());
  pool_type&    ent = conditionPools.back().second;
  iov_key_type& key = conditionPools.back().first;
  key.first         = identifier;
  key.second.first  = make_pair(iov.iovType->name,iov.type);
  key.second.second = iov.key();
  ent               = conditions;
  for(auto c : ent) c->addRef();
  return ent.size();
}

/// Add conditions content to the saved. Note, that dependent conditions shall not be saved!
std::size_t ConditionsTreePersistency::add(const std::string& identifier, ConditionsPool& pool)    {
  DurationStamp stamp(this);
  conditionPools.emplace_back(std::pair<iov_key_type, pool_type>());
  pool_type&    ent = conditionPools.back().second;
  iov_key_type& key = conditionPools.back().first;
  const IOV*    iov = pool.iov;
  key.first         = identifier;
  key.second.first  = make_pair(iov->iovType->name,iov->type);
  key.second.second = iov->key();
  pool.select_all(ent);
  for(auto c : ent) c.ptr()->addRef();
  return ent.size();
}

/// Add conditions content to the saved. Note, that dependent conditions shall not be saved!
std::size_t ConditionsTreePersistency::add(const std::string& identifier, const ConditionsIOVPool& pool)    {
  std::size_t count = 0;
  DurationStamp stamp(this);
  for( const auto& p : pool.elements )  {
    iovPools.emplace_back(std::pair<iov_key_type, pool_type>());
    pool_type&    ent = iovPools.back().second;
    iov_key_type& key = iovPools.back().first;
    const IOV*    iov = p.second->iov;
    key.first         = identifier;
    key.second.first  = make_pair(iov->iovType->name,iov->type);
    key.second.second = p.first;
    p.second->select_all(ent);
    for(auto c : ent) c.ptr()->addRef();
    count += ent.size();
  }
  return count;
}

/// Open ROOT file in read mode
TFile* ConditionsTreePersistency::openFile(const std::string& fname)     {
  TDirectory::TContext context;
  TFile* file = TFile::Open(fname.c_str());
  if ( file && !file->IsZombie()) return file;
  except("ConditionsTreePersistency","+++ FAILED to open ROOT file %s in read-mode.",fname.c_str());
  return 0;
}

/// Clear object content and release allocated memory
void ConditionsTreePersistency::_clear(persistent_type& pool)  {
  /// Cleanup all the stuff not useful....
  for (auto& p : pool )  {
    for(Condition c : p.second )
      c.ptr()->release();
    p.second.clear();
  }
  pool.clear();
}

/// Clear object content and release allocated memory
void ConditionsTreePersistency::clear()  {
  /// Cleanup all the stuff not useful....
  _clear(conditionPools);
  _clear(iovPools);
}

/// Add conditions content to the saved. Note, that dependent conditions shall not be saved!
std::unique_ptr<ConditionsTreePersistency>
ConditionsTreePersistency::load(TFile* file,const std::string& obj)   {
  std::unique_ptr<ConditionsTreePersistency> p;
  if ( file && !file->IsZombie())    {
    TTimeStamp start;
    p.reset((ConditionsTreePersistency*)file->Get(obj.c_str()));
    TTimeStamp stop;
    if ( p.get() )    {
      p->duration = stop.AsDouble()-start.AsDouble();
      return p;
    }
    except("ConditionsTreePersistency",
           "+++ FAILED to load object %s from file %s",
           obj.c_str(), file->GetName());
  }
  except("ConditionsTreePersistency",
         "+++ FAILED to load object %s from file [Invalid file]",obj.c_str());
  return p;
}

/// Load ConditionsPool(s) and populate conditions manager
std::size_t ConditionsTreePersistency::_import(ImportStrategy     strategy,
                                               persistent_type&   persistent_pools,
                                               const std::string& id,
                                               const std::string& iov_type,
                                               const IOV::Key&    iov_key,
                                               ConditionsManager  mgr)
{
  std::size_t count = 0;
  std::pair<bool,const IOVType*> iovTyp(false,0);
  for (auto& iovp : persistent_pools )   {
    bool use = false;
    const iov_key_type& key = iovp.first;
    if ( !(id.empty() || id=="*" || key.first == id) )
      continue;
    if ( !(iov_type.empty() || iov_type == "*" || key.second.first.first == iov_type) )
      continue;
    switch(strategy)   {
    case IMPORT_ALL:
      use = true;
      break;
    case IMPORT_EXACT:
      use = key.second.second == iov_key;
      break;
    case IMPORT_CONTAINED:
      use = IOV::key_is_contained(key.second.second,iov_key);
      break;
    case IMPORT_EDGE_LOWER:
      use = IOV::key_overlaps_lower_end(key.second.second,iov_key);
      break;
    case IMPORT_EDGE_UPPER:
      use = IOV::key_overlaps_higher_end(key.second.second,iov_key);
      break;
    case IMPORT_CONTAINED_LOWER:
      use = IOV::key_is_contained(key.second.second,iov_key) ||
        IOV::key_overlaps_lower_end(key.second.second,iov_key);
      break;
    case IMPORT_CONTAINED_UPPER:
      use = IOV::key_is_contained(key.second.second,iov_key) ||
        IOV::key_overlaps_higher_end(key.second.second,iov_key);
      break;
    default:
      use = false;
      break;
    }
    if ( !use )
      continue;

    iovTyp = mgr.registerIOVType(key.second.first.second,key.second.first.first);
    if ( iovTyp.second )   {
      ConditionsPool* pool = mgr.registerIOV(*iovTyp.second, key.second.second);
      for (Condition c : iovp.second)   {
        Condition::Object* o = c.ptr();
        o->iov = pool->iov;
        if ( pool->insert(o->addRef()) )    {
          //printout(WARNING,"ConditionsTreePersistency","%p: [%016llX] %s -> %s",
          //         o, o->hash,o->GetName(),o->data.str().c_str());
          ++count;
        }
        else   {
          printout(WARNING,"ConditionsTreePersistency",
                   "+++ Ignore condition %s from %s iov:%s [Already present]",
                   c.name(),id.c_str(), iov_type.c_str());
        }
      }
    }
  }
  return count;
}

/// Load ConditionsIOVPool and populate conditions manager
std::size_t ConditionsTreePersistency::importIOVPool(const std::string& identifier,
                                                     const std::string& iov_type,
                                                     ConditionsManager  mgr)
{
  DurationStamp stamp(this);
  return _import(IMPORT_ALL,iovPools,identifier,iov_type,IOV::Key(),mgr);
}

/// Load ConditionsIOVPool and populate conditions manager
std::size_t ConditionsTreePersistency::importConditionsPool(const std::string& identifier,
                                                            const std::string& iov_type,
                                                            ConditionsManager  mgr)
{
  DurationStamp stamp(this);
  return _import(IMPORT_ALL,conditionPools,identifier,iov_type,IOV::Key(),mgr);
}

/// Load conditions pool and populate conditions manager. Allow tro be selective also for the key
std::size_t ConditionsTreePersistency::importConditionsPool(ImportStrategy     strategy,
                                                            const std::string& identifier,
                                                            const std::string& iov_type,
                                                            const IOV::Key&    iov_key,
                                                            ConditionsManager  mgr)   {
  return _import(strategy,conditionPools,identifier,iov_type,iov_key,mgr);
}

/// Save the data content to a root file
int ConditionsTreePersistency::save(TFile* file)    {
  DurationStamp stamp(this);
  //TDirectory::TContext context;
  int nBytes = file->WriteTObject(this,GetName());
  return nBytes;
}

/// Save the data content to a root file
int ConditionsTreePersistency::save(const std::string& fname)    {
  DurationStamp stamp(this);
  //TDirectory::TContext context;
  TFile* file = TFile::Open(fname.c_str(),"RECREATE");
  if ( file && !file->IsZombie())   {
    int nBytes = save(file);
    file->Close();
    delete file;
    return nBytes;
  }
  return -1;
}
