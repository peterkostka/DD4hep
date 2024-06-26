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
#include <DD4hep/OpaqueDataBinder.h>
#include <DD4hep/Conditions.h>
#include <DD4hep/detail/ConditionsInterna.h>

// C/C++ include files
#include <set>
#include <map>
#include <list>
#include <vector>

namespace {

  /// Helper class to bind string values to C++ data objects (primitive or complex)
  template <typename T, typename Q> bool __bind__(const dd4hep::detail::ValueBinder&, T& object, const std::string& val, const Q*)
  {  object.template bind<Q>(val);                 return true;  }

  /// Helper class to bind string values to a STL vector of data objects (primitive or complex)
  template <typename T, typename Q> bool __bind__(const dd4hep::detail::VectorBinder&, T& object, const std::string& val, const Q*)
  {  object.template bind<std::vector<Q> >(val);   return true;  }

  /// Helper class to bind string values to a STL list of data objects (primitive or complex)
  template <typename T, typename Q> bool __bind__(const dd4hep::detail::ListBinder&, T& object, const std::string& val, const Q*)
  {  object.template bind<std::list<Q> >(val);     return true;  }

  /// Helper class to bind string values to a STL set of data objects (primitive or complex)
  template <typename T, typename Q> bool __bind__(const dd4hep::detail::SetBinder&, T& object, const std::string& val, const Q*)
  {  object.template bind<std::set<Q> >(val);      return true;  }

  /// Helper class to bind STL map objects
  template <typename T, typename Q> bool __bind__(const dd4hep::detail::MapBinder&, T& object, const Q*)
  {  object.template bind<Q>();                    return true;  }

}

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {

  /// DD4hep internal namespace declaration for utilities and implementation details
  namespace detail  {

    /// Binding function for scalar items. See the implementation function for the concrete instantiations
    template <typename BINDER, typename T> 
    bool OpaqueDataBinder::bind(const BINDER& b, T& object, const std::string& typ, const std::string& val)  {
#if defined(DD4HEP_HAVE_ALL_PARSERS)
      if ( typ.substr(0,4) == "char" )
        return __bind__(b,object,val,Primitive<char>::null_pointer());
      else if ( typ.substr(0,13) == "unsigned char" )
        return __bind__(b,object,val,Primitive<unsigned char>::null_pointer());
      else if ( typ.substr(0,5) == "short" )
        return __bind__(b,object,val,Primitive<short>::null_pointer());
      else if ( typ.substr(0,14) == "unsigned short" )
        return __bind__(b,object,val,Primitive<unsigned short>::null_pointer());
      else if ( typ.substr(0,12) == "unsigned int" )
        return __bind__(b,object,val,Primitive<unsigned int>::null_pointer());
      else if ( typ.substr(0,13) == "unsigned long" )
        return __bind__(b,object,val,Primitive<unsigned long>::null_pointer());
#else
      // Short and char is not part of the standard dictionaries. Fall back to 'int'.
      if ( typ.substr(0,4) == "char" )
        return __bind__(b,object,val,Primitive<int>::null_pointer());
      else if ( typ.substr(0,5) == "short" )
        return __bind__(b,object,val,Primitive<int>::null_pointer());
#endif
      else if ( typ.substr(0,3) == "int" )
        return __bind__(b,object,val,Primitive<int>::null_pointer());
      else if ( typ.substr(0,4) == "long" ) 
        return __bind__(b,object,val,Primitive<long>::null_pointer());
      else if ( typ.substr(0,5) == "float" )
        return __bind__(b,object,val,Primitive<float>::null_pointer());
      else if ( typ.substr(0,6) == "double" )
        return __bind__(b,object,val,Primitive<double>::null_pointer());
      else if ( typ.substr(0,6) == "string" )
        return __bind__(b,object,val,Primitive<std::string>::null_pointer());
      else if ( typ == "std::string" )
        return __bind__(b,object,val,Primitive<std::string>::null_pointer());
      else if ( typ == "Histo1D" )
        return __bind__(b,object,val,Primitive<std::string>::null_pointer());
      else if ( typ == "Histo2D" )
        return __bind__(b,object,val,Primitive<std::string>::null_pointer());
      else
        printout(INFO,"OpaqueDataBinder","++ Unknown conditions parameter type:%s val:%s",typ.c_str(),val.c_str());
      return __bind__(b,object,val,Primitive<std::string>::null_pointer());
    }
    template bool OpaqueDataBinder::bind<ValueBinder,OpaqueDataBlock>(        const ValueBinder& b, OpaqueDataBlock& object,
                                                                              const std::string& typ, const std::string& val);
    template bool OpaqueDataBinder::bind<VectorBinder,OpaqueDataBlock>(       const VectorBinder& b, OpaqueDataBlock& object,
                                                                              const std::string& typ, const std::string& val);
    template bool OpaqueDataBinder::bind<ListBinder,OpaqueDataBlock>(         const ListBinder& b, OpaqueDataBlock& object,
                                                                              const std::string& typ, const std::string& val);
    template bool OpaqueDataBinder::bind<SetBinder,OpaqueDataBlock>(          const SetBinder& b, OpaqueDataBlock& object,
                                                                              const std::string& typ, const std::string& val);

    template bool OpaqueDataBinder::bind<ValueBinder,Condition>(  const ValueBinder& b, Condition& object,
                                                                  const std::string& typ, const std::string& val);
    template bool OpaqueDataBinder::bind<VectorBinder,Condition>( const VectorBinder& b, Condition& object,
                                                                  const std::string& typ, const std::string& val);
    template bool OpaqueDataBinder::bind<ListBinder,Condition>(   const ListBinder& b, Condition& object,
                                                                  const std::string& typ, const std::string& val);
    template bool OpaqueDataBinder::bind<SetBinder,Condition>(    const SetBinder& b, Condition& object,
                                                                  const std::string& typ, const std::string& val);
  
    /// Binding function for sequences (unmapped STL containers)
    template <typename T> 
    bool OpaqueDataBinder::bind_sequence(T& object, const std::string& typ, const std::string& val)
    {
      std::size_t idx = typ.find('[');
      std::size_t idq = typ.find(']');
      std::string value_type = typ.substr(idx+1,idq-idx-1);
      if ( typ.substr(0,6) == "vector" )
        return bind(VectorBinder(), object, value_type, val);
      else if ( typ.substr(0,4) == "list" )
        return bind(ListBinder(), object, value_type, val);
      else if ( typ.substr(0,3) == "set" )
        return bind(SetBinder(), object, value_type, val);
      else if ( idx == std::string::npos && idq == std::string::npos )
        return bind(ValueBinder(), object, value_type, val);
      return false;
    }
  
    template<typename BINDER, typename OBJECT, typename KEY, typename VAL>
    static void emplace_map_item(const BINDER&,
                                 OBJECT& object,
                                 const KEY& k,
                                 const std::string& val,
                                 const VAL*)
    {
      typedef std::map<KEY,VAL> map_t;
      map_t& m = object.template get<map_t>();
      VAL v;
      if ( !BasicGrammar::instance<VAL>().fromString(&v, val) )  {
        except("OpaqueDataBinder","++ Failed to convert conditions map entry.");
      }
      m.emplace(k,v);
    }

    template<typename BINDER, typename OBJECT, typename KEY> 
    static void emplace_map_key(const BINDER& b,
                                OBJECT& object,
                                const std::string& key_val,
                                const std::string& val_type,
                                const std::string& val,
                                const KEY*)
    {
      KEY key;
      BasicGrammar::instance<KEY>().fromString(&key, key_val);
      // Short and char is not part of the standard dictionaries. Fall back to 'int'.
      if ( val_type.substr(0,4) == "char" )
        emplace_map_item(b, object, key, val, (int*)0);
      else if ( val_type.substr(0,5) == "short" )
        emplace_map_item(b, object, key, val, (int*)0);
      else if ( val_type.substr(0,3) == "int" )
        emplace_map_item(b, object, key, val, (int*)0);
      else if ( val_type.substr(0,4) == "long" )
        emplace_map_item(b, object, key, val, (long*)0);
      else if ( val_type.substr(0,5) == "float" )
        emplace_map_item(b, object, key, val, (float*)0);
      else if ( val_type.substr(0,6) == "double" )
        emplace_map_item(b, object, key, val, (double*)0);
      else if ( val_type.substr(0,6) == "string" )
        emplace_map_item(b, object, key, val, (std::string*)0);
      else if ( val_type == "std::string" )
        emplace_map_item(b, object, key, val, (std::string*)0);
      else {
        printout(INFO,"Param","++ Unknown conditions parameter type:%s data:%s",
                 val_type.c_str(),val.c_str());
        emplace_map_item(b, object, key, val, (std::string*)0);
      }
    }

    template<typename BINDER, typename OBJECT, typename KEY, typename VAL>
    static void emplace_map_pair(const BINDER&,
                                 OBJECT& object,
                                 const std::string& data,
                                 const KEY*,
                                 const VAL*)
    {
      typedef std::map<KEY,VAL> map_t;
      std::pair<KEY,VAL> entry;
      map_t& m = object.template get<map_t>();
      if ( !BasicGrammar::instance<std::pair<KEY,VAL> >().fromString(&entry,data) )  {
        except("OpaqueDataBinder","++ Failed to convert conditions map entry.");
      }
      m.emplace(entry);
    }

    template<typename BINDER, typename OBJECT, typename KEY> 
    static void emplace_map_data(const BINDER& b,
                                 OBJECT& object,
                                 const std::string& val_type,
                                 const std::string& pair_data,
                                 const KEY*)
    {
      // Short and char is not part of the standard dictionaries. Fall back to 'int'.
      if ( val_type.substr(0,4) == "char" )
        emplace_map_pair(b, object, pair_data, (KEY*)0, (int*)0);
      else if ( val_type.substr(0,5) == "short" )
        emplace_map_pair(b, object, pair_data, (KEY*)0, (int*)0);
      else if ( val_type.substr(0,3) == "int" )
        emplace_map_pair(b, object, pair_data, (KEY*)0, (int*)0);
      else if ( val_type.substr(0,4) == "long" )
        emplace_map_pair(b, object, pair_data, (KEY*)0, (long*)0);
      else if ( val_type.substr(0,5) == "float" )
        emplace_map_pair(b, object, pair_data, (KEY*)0, (float*)0);
      else if ( val_type.substr(0,6) == "double" )
        emplace_map_pair(b, object, pair_data, (KEY*)0, (double*)0);
      else if ( val_type.substr(0,6) == "string" )
        emplace_map_pair(b, object, pair_data, (KEY*)0, (std::string*)0);
      else if ( val_type == "std::string" )
        emplace_map_pair(b, object, pair_data, (KEY*)0, (std::string*)0);
      else {
        printout(INFO,"Param","++ Unknown conditions parameter type:%s data:%s",
                 val_type.c_str(),pair_data.c_str());
        emplace_map_pair(b, object, pair_data, (KEY*)0, (std::string*)0);
      }
    }

    template<typename BINDER, typename OBJECT, typename KEY> 
    static void bind_mapping(const BINDER& b, const std::string& val_type, OBJECT& object, const KEY*)   {
      if ( val_type.substr(0,3) == "int" )
        __bind__(b,object, (std::map<KEY,int>*)0);
#if defined(DD4HEP_HAVE_ALL_PARSERS)
      else if ( val_type.substr(0,12) == "unsigned int" )
        __bind__(b,object, (std::map<KEY,unsigned int>*)0);
      else if ( val_type.substr(0,4) == "char" )
        __bind__(b,object, (std::map<KEY,char>*)0);
      else if ( val_type.substr(0,13) == "unsigned char" )
        __bind__(b,object, (std::map<KEY,unsigned char>*)0);
      else if ( val_type.substr(0,5) == "short" )
        __bind__(b,object, (std::map<KEY,short>*)0);
      else if ( val_type.substr(0,14) == "unsigned short" )
        __bind__(b,object, (std::map<KEY,unsigned short>*)0);
      else if ( val_type.substr(0,13) == "unsigned long" )
        __bind__(b,object, (std::map<KEY,unsigned long>*)0);
#else
      // Short and char is not part of the standard dictionaries. Fall back to 'int'.
      else if ( val_type.substr(0,4) == "char" )
        __bind__(b,object, (std::map<KEY,int>*)0);
      else if ( val_type.substr(0,5) == "short" )
        __bind__(b,object, (std::map<KEY,int>*)0);
#endif
      else if ( val_type.substr(0,4) == "long" )
        __bind__(b,object, (std::map<KEY,long>*)0);
      else if ( val_type.substr(0,5) == "float" )
        __bind__(b,object, (std::map<KEY,float>*)0);
      else if ( val_type.substr(0,6) == "double" )
        __bind__(b,object, (std::map<KEY,double>*)0);
      else if ( val_type.substr(0,6) == "string" )
        __bind__(b,object, (std::map<KEY,std::string>*)0);
      else if ( val_type == "std::string" )
        __bind__(b,object, (std::map<KEY,std::string>*)0);
      else {
        __bind__(b,object, (std::map<KEY,std::string>*)0);
      }
    }
  
    /// Binding function for STL maps
    template <typename BINDER, typename OBJECT> 
    bool OpaqueDataBinder::bind_map(const BINDER& b, OBJECT& object,
                                    const std::string& key_type, const std::string& val_type)   {
      // Short and char is not part of the standard dictionaries. Fall back to 'int'.
      if ( key_type.substr(0,3) == "int" )
        bind_mapping(b, val_type, object, Primitive<int>::null_pointer());
#if defined(DD4HEP_HAVE_ALL_PARSERS)
      else if ( key_type.substr(0,4) == "char" )
        bind_mapping(b, val_type, object, Primitive<int>::null_pointer());
      else if ( key_type.substr(0,5) == "short" )
        bind_mapping(b, val_type, object, Primitive<int>::null_pointer());
      else if ( key_type.substr(0,4) == "long" )
        bind_mapping(b, val_type, object, Primitive<long>::null_pointer());
      else if ( key_type.substr(0,5) == "float" )
        bind_mapping(b, val_type, object, Primitive<float>::null_pointer());
      else if ( key_type.substr(0,6) == "double" )
        bind_mapping(b, val_type, object, Primitive<double>::null_pointer());
#endif
      else if ( key_type.substr(0,6) == "string" )
        bind_mapping(b, val_type, object, Primitive<std::string>::null_pointer());
      else if ( key_type == "std::string" )
        bind_mapping(b, val_type, object, Primitive<std::string>::null_pointer());
      else {
        printout(INFO,"OpaqueDataBinder","++ Unknown MAP-conditions key-type:%s",key_type.c_str());
        bind_mapping(b, val_type, object, Primitive<std::string>::null_pointer());
      }
      return true;
    }

    /// Filling function for STL maps.
    template <typename BINDER, typename OBJECT>
    bool OpaqueDataBinder::insert_map(const BINDER& b, OBJECT& object,
                                      const std::string& key_type, const std::string& key,
                                      const std::string& val_type, const std::string& val)
    {
      if ( key_type.substr(0,3) == "int" )
        emplace_map_key(b, object, key, val_type, val, Primitive<int>::null_pointer());
#if defined(DD4HEP_HAVE_ALL_PARSERS)
      // Short and char is not part of the standard dictionaries. Fall back to 'int'.
      else if ( key_type.substr(0,4) == "char" )
        emplace_map_key(b, object, key, val_type, val, Primitive<int>::null_pointer());
      else if ( key_type.substr(0,5) == "short" )
        emplace_map_key(b, object, key, val_type, val, Primitive<int>::null_pointer());
      else if ( key_type.substr(0,4) == "long" )
        emplace_map_key(b, object, key, val_type, val, Primitive<long>::null_pointer());
      else if ( key_type.substr(0,5) == "float" )
        emplace_map_key(b, object, key, val_type, val, Primitive<float>::null_pointer());
      else if ( key_type.substr(0,6) == "double" )
        emplace_map_key(b, object, key, val_type, val, Primitive<double>::null_pointer());
#endif
      else if ( key_type.substr(0,6) == "string" )
        emplace_map_key(b, object, key, val_type, val, Primitive<std::string>::null_pointer());
      else if ( key_type == "std::string" )
        emplace_map_key(b, object, key, val_type, val, Primitive<std::string>::null_pointer());
      else {
        printout(INFO,"OpaqueDataBinder","++ Unknown MAP-conditions key-type:%s",key_type.c_str());
        emplace_map_key(b, object, key, val_type, val, Primitive<std::string>::null_pointer());
      }
      return true;
    }

    /// Filling function for STL maps.
    template <typename BINDER, typename OBJECT> 
    bool OpaqueDataBinder::insert_map(const BINDER& b, OBJECT& object,
                                      const std::string& key_type, const std::string& val_type,
                                      const std::string& pair_data)
    {
      if ( key_type.substr(0,3) == "int" )
        emplace_map_data(b, object, val_type, pair_data, Primitive<int>::null_pointer());
#if defined(DD4HEP_HAVE_ALL_PARSERS)
      // Short and char is not part of the standard dictionaries. Fall back to 'int'.
      else if ( key_type.substr(0,4) == "char" )
        emplace_map_data(b, object, val_type, pair_data, Primitive<int>::null_pointer());
      else if ( key_type.substr(0,5) == "short" )
        emplace_map_data(b, object, val_type, pair_data, Primitive<int>::null_pointer());
      else if ( key_type.substr(0,4) == "long" )
        emplace_map_data(b, object, val_type, pair_data, Primitive<long>::null_pointer());
      else if ( key_type.substr(0,5) == "float" )
        emplace_map_data(b, object, val_type, pair_data, Primitive<float>::null_pointer());
      else if ( key_type.substr(0,6) == "double" )
        emplace_map_data(b, object, val_type, pair_data, Primitive<double>::null_pointer());
#endif
      else if ( key_type.substr(0,6) == "string" )
        emplace_map_data(b, object, val_type, pair_data, Primitive<std::string>::null_pointer());
      else if ( key_type == "std::string" )
        emplace_map_data(b, object, val_type, pair_data, Primitive<std::string>::null_pointer());
      else {
        printout(INFO,"OpaqueDataBinder","++ Unknown MAP-conditions key-type:%s",key_type.c_str());
        emplace_map_data(b, object, val_type, pair_data, Primitive<std::string>::null_pointer());
      }
      return true;
    }

    template bool OpaqueDataBinder::bind_sequence<OpaqueDataBlock>(OpaqueDataBlock& object,const std::string& typ,const std::string& val);
    template bool OpaqueDataBinder::bind_map<MapBinder,OpaqueDataBlock>(    const MapBinder& b, OpaqueDataBlock& object,
                                                                            const std::string& typ,const std::string& val);
    template bool OpaqueDataBinder::insert_map<MapBinder,OpaqueDataBlock>(  const MapBinder& b, OpaqueDataBlock& object,
                                                                            const std::string& key_type, const std::string& key,
                                                                            const std::string& val_type, const std::string& val);
    template bool OpaqueDataBinder::insert_map<MapBinder,OpaqueDataBlock>(  const MapBinder& b, OpaqueDataBlock& object,
                                                                            const std::string& key_type, const std::string& val_type,
                                                                            const std::string& pair_data);

    /// Instantiation for Conditions:
    template bool OpaqueDataBinder::bind_sequence<Condition>(   Condition& object,
                                                                const std::string& typ,const std::string& val);
    /// Conditions binding function for STL maps
    template <> bool OpaqueDataBinder::bind_map(const MapBinder& b, Condition& object,
                                                const std::string& key_type, const std::string& val_type)
    {    return bind_map(b, object->data, key_type, val_type);  }

    /// Conditions: Filling function for STL maps.
    template <> bool OpaqueDataBinder::insert_map(const MapBinder& b, Condition& object,
                                                  const std::string& key_type, const std::string& key,
                                                  const std::string& val_type, const std::string& val)
    {    return insert_map(b, object->data, key_type, key, val_type, val);    }

    /// Conditions: Filling function for STL maps.
    template <> bool OpaqueDataBinder::insert_map(const MapBinder& b, Condition& object,
                                                  const std::string& key_type, const std::string& val_type, const std::string& pair_data)
    {    return insert_map(b, object->data, key_type, val_type, pair_data);   }  
  }
}
