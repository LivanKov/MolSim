// Copyright (c) 2005-2014 Code Synthesis Tools CC
//
// This program was generated by CodeSynthesis XSD, an XML Schema to
// C++ data binding compiler.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 2 as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
//
// In addition, as a special exception, Code Synthesis Tools CC gives
// permission to link this program with the Xerces-C++ library (or with
// modified versions of Xerces-C++ that use the same license as Xerces-C++),
// and distribute linked combinations including the two. You must obey
// the GNU General Public License version 2 in all respects for all of
// the code used other than Xerces-C++. If you modify this copy of the
// program, you may extend this exception to your version of the program,
// but you are not obligated to do so. If you do not wish to do so, delete
// this exception statement from your version.
//
// Furthermore, Code Synthesis Tools CC makes a special exception for
// the Free/Libre and Open Source Software (FLOSS) which is described
// in the accompanying FLOSSE file.
//

#ifndef MOL_SIM_HXX
#define MOL_SIM_HXX

#ifndef XSD_USE_CHAR
#define XSD_USE_CHAR
#endif

#ifndef XSD_CXX_TREE_USE_CHAR
#define XSD_CXX_TREE_USE_CHAR
#endif

// Begin prologue.
//
//
// End prologue.

#include <xsd/cxx/config.hxx>

#if (XSD_INT_VERSION != 4000000L)
#error XSD runtime version mismatch
#endif

#include <xsd/cxx/pre.hxx>

#include <xsd/cxx/xml/char-utf8.hxx>

#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/elements.hxx>
#include <xsd/cxx/tree/types.hxx>

#include <xsd/cxx/xml/error-handler.hxx>

#include <xsd/cxx/xml/dom/auto-ptr.hxx>

#include <xsd/cxx/tree/parsing.hxx>
#include <xsd/cxx/tree/parsing/byte.hxx>
#include <xsd/cxx/tree/parsing/unsigned-byte.hxx>
#include <xsd/cxx/tree/parsing/short.hxx>
#include <xsd/cxx/tree/parsing/unsigned-short.hxx>
#include <xsd/cxx/tree/parsing/int.hxx>
#include <xsd/cxx/tree/parsing/unsigned-int.hxx>
#include <xsd/cxx/tree/parsing/long.hxx>
#include <xsd/cxx/tree/parsing/unsigned-long.hxx>
#include <xsd/cxx/tree/parsing/boolean.hxx>
#include <xsd/cxx/tree/parsing/float.hxx>
#include <xsd/cxx/tree/parsing/double.hxx>
#include <xsd/cxx/tree/parsing/decimal.hxx>

namespace xml_schema
{
  // anyType and anySimpleType.
  //
  typedef ::xsd::cxx::tree::type type;
  typedef ::xsd::cxx::tree::simple_type< char, type > simple_type;
  typedef ::xsd::cxx::tree::type container;

  // 8-bit
  //
  typedef signed char byte;
  typedef unsigned char unsigned_byte;

  // 16-bit
  //
  typedef short short_;
  typedef unsigned short unsigned_short;

  // 32-bit
  //
  typedef int int_;
  typedef unsigned int unsigned_int;

  // 64-bit
  //
  typedef long long long_;
  typedef unsigned long long unsigned_long;

  // Supposed to be arbitrary-length integral types.
  //
  typedef long long integer;
  typedef long long non_positive_integer;
  typedef unsigned long long non_negative_integer;
  typedef unsigned long long positive_integer;
  typedef long long negative_integer;

  // Boolean.
  //
  typedef bool boolean;

  // Floating-point types.
  //
  typedef float float_;
  typedef double double_;
  typedef double decimal;

  // String types.
  //
  typedef ::xsd::cxx::tree::string< char, simple_type > string;
  typedef ::xsd::cxx::tree::normalized_string< char, string > normalized_string;
  typedef ::xsd::cxx::tree::token< char, normalized_string > token;
  typedef ::xsd::cxx::tree::name< char, token > name;
  typedef ::xsd::cxx::tree::nmtoken< char, token > nmtoken;
  typedef ::xsd::cxx::tree::nmtokens< char, simple_type, nmtoken > nmtokens;
  typedef ::xsd::cxx::tree::ncname< char, name > ncname;
  typedef ::xsd::cxx::tree::language< char, token > language;

  // ID/IDREF.
  //
  typedef ::xsd::cxx::tree::id< char, ncname > id;
  typedef ::xsd::cxx::tree::idref< char, ncname, type > idref;
  typedef ::xsd::cxx::tree::idrefs< char, simple_type, idref > idrefs;

  // URI.
  //
  typedef ::xsd::cxx::tree::uri< char, simple_type > uri;

  // Qualified name.
  //
  typedef ::xsd::cxx::tree::qname< char, simple_type, uri, ncname > qname;

  // Binary.
  //
  typedef ::xsd::cxx::tree::buffer< char > buffer;
  typedef ::xsd::cxx::tree::base64_binary< char, simple_type > base64_binary;
  typedef ::xsd::cxx::tree::hex_binary< char, simple_type > hex_binary;

  // Date/time.
  //
  typedef ::xsd::cxx::tree::time_zone time_zone;
  typedef ::xsd::cxx::tree::date< char, simple_type > date;
  typedef ::xsd::cxx::tree::date_time< char, simple_type > date_time;
  typedef ::xsd::cxx::tree::duration< char, simple_type > duration;
  typedef ::xsd::cxx::tree::gday< char, simple_type > gday;
  typedef ::xsd::cxx::tree::gmonth< char, simple_type > gmonth;
  typedef ::xsd::cxx::tree::gmonth_day< char, simple_type > gmonth_day;
  typedef ::xsd::cxx::tree::gyear< char, simple_type > gyear;
  typedef ::xsd::cxx::tree::gyear_month< char, simple_type > gyear_month;
  typedef ::xsd::cxx::tree::time< char, simple_type > time;

  // Entity.
  //
  typedef ::xsd::cxx::tree::entity< char, ncname > entity;
  typedef ::xsd::cxx::tree::entities< char, simple_type, entity > entities;

  typedef ::xsd::cxx::tree::content_order content_order;
  // Flags and properties.
  //
  typedef ::xsd::cxx::tree::flags flags;
  typedef ::xsd::cxx::tree::properties< char > properties;

  // Parsing/serialization diagnostics.
  //
  typedef ::xsd::cxx::tree::severity severity;
  typedef ::xsd::cxx::tree::error< char > error;
  typedef ::xsd::cxx::tree::diagnostics< char > diagnostics;

  // Exceptions.
  //
  typedef ::xsd::cxx::tree::exception< char > exception;
  typedef ::xsd::cxx::tree::bounds< char > bounds;
  typedef ::xsd::cxx::tree::duplicate_id< char > duplicate_id;
  typedef ::xsd::cxx::tree::parsing< char > parsing;
  typedef ::xsd::cxx::tree::expected_element< char > expected_element;
  typedef ::xsd::cxx::tree::unexpected_element< char > unexpected_element;
  typedef ::xsd::cxx::tree::expected_attribute< char > expected_attribute;
  typedef ::xsd::cxx::tree::unexpected_enumerator< char > unexpected_enumerator;
  typedef ::xsd::cxx::tree::expected_text_content< char > expected_text_content;
  typedef ::xsd::cxx::tree::no_prefix_mapping< char > no_prefix_mapping;

  // Error handler callback interface.
  //
  typedef ::xsd::cxx::xml::error_handler< char > error_handler;

  // DOM interaction.
  //
  namespace dom
  {
    // Automatic pointer for DOMDocument.
    //
    using ::xsd::cxx::xml::dom::auto_ptr;

#ifndef XSD_CXX_TREE_TREE_NODE_KEY__XML_SCHEMA
#define XSD_CXX_TREE_TREE_NODE_KEY__XML_SCHEMA
    // DOM user data key for back pointers to tree nodes.
    //
    const XMLCh* const tree_node_key = ::xsd::cxx::tree::user_data_keys::node;
#endif
  }
}

// Forward declarations.
//
class MolSim;
class simulation_parameters;
class discs;
class cuboids;
class domain_size;
class disc;
class cuboid;
class center;
class initial_velocity;
class coordinate;
class dimensions;

#include <memory>    // ::std::auto_ptr
#include <limits>    // std::numeric_limits
#include <algorithm> // std::binary_search

#include <xsd/cxx/xml/char-utf8.hxx>

#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/elements.hxx>
#include <xsd/cxx/tree/containers.hxx>
#include <xsd/cxx/tree/list.hxx>

#include <xsd/cxx/xml/dom/parsing-header.hxx>

class MolSim: public ::xml_schema::type
{
  public:
  // simulation_parameters
  //
  typedef ::simulation_parameters simulation_parameters_type;
  typedef ::xsd::cxx::tree::traits< simulation_parameters_type, char > simulation_parameters_traits;

  const simulation_parameters_type&
  simulation_parameters () const;

  simulation_parameters_type&
  simulation_parameters ();

  void
  simulation_parameters (const simulation_parameters_type& x);

  void
  simulation_parameters (::std::auto_ptr< simulation_parameters_type > p);

  // discs
  //
  typedef ::discs discs_type;
  typedef ::xsd::cxx::tree::traits< discs_type, char > discs_traits;

  const discs_type&
  discs () const;

  discs_type&
  discs ();

  void
  discs (const discs_type& x);

  void
  discs (::std::auto_ptr< discs_type > p);

  // cuboids
  //
  typedef ::cuboids cuboids_type;
  typedef ::xsd::cxx::tree::traits< cuboids_type, char > cuboids_traits;

  const cuboids_type&
  cuboids () const;

  cuboids_type&
  cuboids ();

  void
  cuboids (const cuboids_type& x);

  void
  cuboids (::std::auto_ptr< cuboids_type > p);

  // Constructors.
  //
  MolSim (const simulation_parameters_type&,
          const discs_type&,
          const cuboids_type&);

  MolSim (::std::auto_ptr< simulation_parameters_type >,
          ::std::auto_ptr< discs_type >,
          ::std::auto_ptr< cuboids_type >);

  MolSim (const ::xercesc::DOMElement& e,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  MolSim (const MolSim& x,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  virtual MolSim*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  MolSim&
  operator= (const MolSim& x);

  virtual 
  ~MolSim ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< simulation_parameters_type > simulation_parameters_;
  ::xsd::cxx::tree::one< discs_type > discs_;
  ::xsd::cxx::tree::one< cuboids_type > cuboids_;
};

class simulation_parameters: public ::xml_schema::type
{
  public:
  // end_time
  //
  typedef ::xml_schema::double_ end_time_type;
  typedef ::xsd::cxx::tree::traits< end_time_type, char, ::xsd::cxx::tree::schema_type::double_ > end_time_traits;

  const end_time_type&
  end_time () const;

  end_time_type&
  end_time ();

  void
  end_time (const end_time_type& x);

  // delta_time
  //
  typedef ::xml_schema::double_ delta_time_type;
  typedef ::xsd::cxx::tree::traits< delta_time_type, char, ::xsd::cxx::tree::schema_type::double_ > delta_time_traits;

  const delta_time_type&
  delta_time () const;

  delta_time_type&
  delta_time ();

  void
  delta_time (const delta_time_type& x);

  // output_basename
  //
  typedef ::xml_schema::string output_basename_type;
  typedef ::xsd::cxx::tree::traits< output_basename_type, char > output_basename_traits;

  const output_basename_type&
  output_basename () const;

  output_basename_type&
  output_basename ();

  void
  output_basename (const output_basename_type& x);

  void
  output_basename (::std::auto_ptr< output_basename_type > p);

  // write_frequency
  //
  typedef ::xml_schema::unsigned_int write_frequency_type;
  typedef ::xsd::cxx::tree::traits< write_frequency_type, char > write_frequency_traits;

  const write_frequency_type&
  write_frequency () const;

  write_frequency_type&
  write_frequency ();

  void
  write_frequency (const write_frequency_type& x);

  // r_cutoff_radius
  //
  typedef ::xml_schema::double_ r_cutoff_radius_type;
  typedef ::xsd::cxx::tree::traits< r_cutoff_radius_type, char, ::xsd::cxx::tree::schema_type::double_ > r_cutoff_radius_traits;

  const r_cutoff_radius_type&
  r_cutoff_radius () const;

  r_cutoff_radius_type&
  r_cutoff_radius ();

  void
  r_cutoff_radius (const r_cutoff_radius_type& x);

  // domain_size
  //
  typedef ::domain_size domain_size_type;
  typedef ::xsd::cxx::tree::traits< domain_size_type, char > domain_size_traits;

  const domain_size_type&
  domain_size () const;

  domain_size_type&
  domain_size ();

  void
  domain_size (const domain_size_type& x);

  void
  domain_size (::std::auto_ptr< domain_size_type > p);

  // Constructors.
  //
  simulation_parameters (const end_time_type&,
                         const delta_time_type&,
                         const output_basename_type&,
                         const write_frequency_type&,
                         const r_cutoff_radius_type&,
                         const domain_size_type&);

  simulation_parameters (const end_time_type&,
                         const delta_time_type&,
                         const output_basename_type&,
                         const write_frequency_type&,
                         const r_cutoff_radius_type&,
                         ::std::auto_ptr< domain_size_type >);

  simulation_parameters (const ::xercesc::DOMElement& e,
                         ::xml_schema::flags f = 0,
                         ::xml_schema::container* c = 0);

  simulation_parameters (const simulation_parameters& x,
                         ::xml_schema::flags f = 0,
                         ::xml_schema::container* c = 0);

  virtual simulation_parameters*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  simulation_parameters&
  operator= (const simulation_parameters& x);

  virtual 
  ~simulation_parameters ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< end_time_type > end_time_;
  ::xsd::cxx::tree::one< delta_time_type > delta_time_;
  ::xsd::cxx::tree::one< output_basename_type > output_basename_;
  ::xsd::cxx::tree::one< write_frequency_type > write_frequency_;
  ::xsd::cxx::tree::one< r_cutoff_radius_type > r_cutoff_radius_;
  ::xsd::cxx::tree::one< domain_size_type > domain_size_;
};

class discs: public ::xml_schema::type
{
  public:
  // disc
  //
  typedef ::disc disc_type;
  typedef ::xsd::cxx::tree::sequence< disc_type > disc_sequence;
  typedef disc_sequence::iterator disc_iterator;
  typedef disc_sequence::const_iterator disc_const_iterator;
  typedef ::xsd::cxx::tree::traits< disc_type, char > disc_traits;

  const disc_sequence&
  disc () const;

  disc_sequence&
  disc ();

  void
  disc (const disc_sequence& s);

  // Constructors.
  //
  discs ();

  discs (const ::xercesc::DOMElement& e,
         ::xml_schema::flags f = 0,
         ::xml_schema::container* c = 0);

  discs (const discs& x,
         ::xml_schema::flags f = 0,
         ::xml_schema::container* c = 0);

  virtual discs*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  discs&
  operator= (const discs& x);

  virtual 
  ~discs ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  disc_sequence disc_;
};

class cuboids: public ::xml_schema::type
{
  public:
  // cuboid
  //
  typedef ::cuboid cuboid_type;
  typedef ::xsd::cxx::tree::sequence< cuboid_type > cuboid_sequence;
  typedef cuboid_sequence::iterator cuboid_iterator;
  typedef cuboid_sequence::const_iterator cuboid_const_iterator;
  typedef ::xsd::cxx::tree::traits< cuboid_type, char > cuboid_traits;

  const cuboid_sequence&
  cuboid () const;

  cuboid_sequence&
  cuboid ();

  void
  cuboid (const cuboid_sequence& s);

  // Constructors.
  //
  cuboids ();

  cuboids (const ::xercesc::DOMElement& e,
           ::xml_schema::flags f = 0,
           ::xml_schema::container* c = 0);

  cuboids (const cuboids& x,
           ::xml_schema::flags f = 0,
           ::xml_schema::container* c = 0);

  virtual cuboids*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  cuboids&
  operator= (const cuboids& x);

  virtual 
  ~cuboids ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  cuboid_sequence cuboid_;
};

class domain_size: public ::xml_schema::type
{
  public:
  // x
  //
  typedef ::xml_schema::double_ x_type;
  typedef ::xsd::cxx::tree::traits< x_type, char, ::xsd::cxx::tree::schema_type::double_ > x_traits;

  const x_type&
  x () const;

  x_type&
  x ();

  void
  x (const x_type& x);

  // y
  //
  typedef ::xml_schema::double_ y_type;
  typedef ::xsd::cxx::tree::traits< y_type, char, ::xsd::cxx::tree::schema_type::double_ > y_traits;

  const y_type&
  y () const;

  y_type&
  y ();

  void
  y (const y_type& x);

  // z
  //
  typedef ::xml_schema::double_ z_type;
  typedef ::xsd::cxx::tree::traits< z_type, char, ::xsd::cxx::tree::schema_type::double_ > z_traits;

  const z_type&
  z () const;

  z_type&
  z ();

  void
  z (const z_type& x);

  // Constructors.
  //
  domain_size (const x_type&,
               const y_type&,
               const z_type&);

  domain_size (const ::xercesc::DOMElement& e,
               ::xml_schema::flags f = 0,
               ::xml_schema::container* c = 0);

  domain_size (const domain_size& x,
               ::xml_schema::flags f = 0,
               ::xml_schema::container* c = 0);

  virtual domain_size*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  domain_size&
  operator= (const domain_size& x);

  virtual 
  ~domain_size ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< x_type > x_;
  ::xsd::cxx::tree::one< y_type > y_;
  ::xsd::cxx::tree::one< z_type > z_;
};

class disc: public ::xml_schema::type
{
  public:
  // center
  //
  typedef ::center center_type;
  typedef ::xsd::cxx::tree::traits< center_type, char > center_traits;

  const center_type&
  center () const;

  center_type&
  center ();

  void
  center (const center_type& x);

  void
  center (::std::auto_ptr< center_type > p);

  // initial_velocity
  //
  typedef ::initial_velocity initial_velocity_type;
  typedef ::xsd::cxx::tree::traits< initial_velocity_type, char > initial_velocity_traits;

  const initial_velocity_type&
  initial_velocity () const;

  initial_velocity_type&
  initial_velocity ();

  void
  initial_velocity (const initial_velocity_type& x);

  void
  initial_velocity (::std::auto_ptr< initial_velocity_type > p);

  // radius
  //
  typedef ::xml_schema::unsigned_int radius_type;
  typedef ::xsd::cxx::tree::traits< radius_type, char > radius_traits;

  const radius_type&
  radius () const;

  radius_type&
  radius ();

  void
  radius (const radius_type& x);

  // mesh_width
  //
  typedef ::xml_schema::double_ mesh_width_type;
  typedef ::xsd::cxx::tree::traits< mesh_width_type, char, ::xsd::cxx::tree::schema_type::double_ > mesh_width_traits;

  const mesh_width_type&
  mesh_width () const;

  mesh_width_type&
  mesh_width ();

  void
  mesh_width (const mesh_width_type& x);

  // mass
  //
  typedef ::xml_schema::double_ mass_type;
  typedef ::xsd::cxx::tree::traits< mass_type, char, ::xsd::cxx::tree::schema_type::double_ > mass_traits;

  const mass_type&
  mass () const;

  mass_type&
  mass ();

  void
  mass (const mass_type& x);

  // Constructors.
  //
  disc (const center_type&,
        const initial_velocity_type&,
        const radius_type&,
        const mesh_width_type&,
        const mass_type&);

  disc (::std::auto_ptr< center_type >,
        ::std::auto_ptr< initial_velocity_type >,
        const radius_type&,
        const mesh_width_type&,
        const mass_type&);

  disc (const ::xercesc::DOMElement& e,
        ::xml_schema::flags f = 0,
        ::xml_schema::container* c = 0);

  disc (const disc& x,
        ::xml_schema::flags f = 0,
        ::xml_schema::container* c = 0);

  virtual disc*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  disc&
  operator= (const disc& x);

  virtual 
  ~disc ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< center_type > center_;
  ::xsd::cxx::tree::one< initial_velocity_type > initial_velocity_;
  ::xsd::cxx::tree::one< radius_type > radius_;
  ::xsd::cxx::tree::one< mesh_width_type > mesh_width_;
  ::xsd::cxx::tree::one< mass_type > mass_;
};

class cuboid: public ::xml_schema::type
{
  public:
  // coordinate
  //
  typedef ::coordinate coordinate_type;
  typedef ::xsd::cxx::tree::traits< coordinate_type, char > coordinate_traits;

  const coordinate_type&
  coordinate () const;

  coordinate_type&
  coordinate ();

  void
  coordinate (const coordinate_type& x);

  void
  coordinate (::std::auto_ptr< coordinate_type > p);

  // dimensions
  //
  typedef ::dimensions dimensions_type;
  typedef ::xsd::cxx::tree::traits< dimensions_type, char > dimensions_traits;

  const dimensions_type&
  dimensions () const;

  dimensions_type&
  dimensions ();

  void
  dimensions (const dimensions_type& x);

  void
  dimensions (::std::auto_ptr< dimensions_type > p);

  // mesh_width
  //
  typedef ::xml_schema::double_ mesh_width_type;
  typedef ::xsd::cxx::tree::traits< mesh_width_type, char, ::xsd::cxx::tree::schema_type::double_ > mesh_width_traits;

  const mesh_width_type&
  mesh_width () const;

  mesh_width_type&
  mesh_width ();

  void
  mesh_width (const mesh_width_type& x);

  // mass
  //
  typedef ::xml_schema::double_ mass_type;
  typedef ::xsd::cxx::tree::traits< mass_type, char, ::xsd::cxx::tree::schema_type::double_ > mass_traits;

  const mass_type&
  mass () const;

  mass_type&
  mass ();

  void
  mass (const mass_type& x);

  // initial_velocity
  //
  typedef ::initial_velocity initial_velocity_type;
  typedef ::xsd::cxx::tree::traits< initial_velocity_type, char > initial_velocity_traits;

  const initial_velocity_type&
  initial_velocity () const;

  initial_velocity_type&
  initial_velocity ();

  void
  initial_velocity (const initial_velocity_type& x);

  void
  initial_velocity (::std::auto_ptr< initial_velocity_type > p);

  // average_velocity
  //
  typedef ::xml_schema::double_ average_velocity_type;
  typedef ::xsd::cxx::tree::traits< average_velocity_type, char, ::xsd::cxx::tree::schema_type::double_ > average_velocity_traits;

  const average_velocity_type&
  average_velocity () const;

  average_velocity_type&
  average_velocity ();

  void
  average_velocity (const average_velocity_type& x);

  // Constructors.
  //
  cuboid (const coordinate_type&,
          const dimensions_type&,
          const mesh_width_type&,
          const mass_type&,
          const initial_velocity_type&,
          const average_velocity_type&);

  cuboid (::std::auto_ptr< coordinate_type >,
          ::std::auto_ptr< dimensions_type >,
          const mesh_width_type&,
          const mass_type&,
          ::std::auto_ptr< initial_velocity_type >,
          const average_velocity_type&);

  cuboid (const ::xercesc::DOMElement& e,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  cuboid (const cuboid& x,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  virtual cuboid*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  cuboid&
  operator= (const cuboid& x);

  virtual 
  ~cuboid ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< coordinate_type > coordinate_;
  ::xsd::cxx::tree::one< dimensions_type > dimensions_;
  ::xsd::cxx::tree::one< mesh_width_type > mesh_width_;
  ::xsd::cxx::tree::one< mass_type > mass_;
  ::xsd::cxx::tree::one< initial_velocity_type > initial_velocity_;
  ::xsd::cxx::tree::one< average_velocity_type > average_velocity_;
};

class center: public ::xml_schema::type
{
  public:
  // x
  //
  typedef ::xml_schema::double_ x_type;
  typedef ::xsd::cxx::tree::traits< x_type, char, ::xsd::cxx::tree::schema_type::double_ > x_traits;

  const x_type&
  x () const;

  x_type&
  x ();

  void
  x (const x_type& x);

  // y
  //
  typedef ::xml_schema::double_ y_type;
  typedef ::xsd::cxx::tree::traits< y_type, char, ::xsd::cxx::tree::schema_type::double_ > y_traits;

  const y_type&
  y () const;

  y_type&
  y ();

  void
  y (const y_type& x);

  // z
  //
  typedef ::xml_schema::double_ z_type;
  typedef ::xsd::cxx::tree::traits< z_type, char, ::xsd::cxx::tree::schema_type::double_ > z_traits;

  const z_type&
  z () const;

  z_type&
  z ();

  void
  z (const z_type& x);

  // Constructors.
  //
  center (const x_type&,
          const y_type&,
          const z_type&);

  center (const ::xercesc::DOMElement& e,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  center (const center& x,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  virtual center*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  center&
  operator= (const center& x);

  virtual 
  ~center ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< x_type > x_;
  ::xsd::cxx::tree::one< y_type > y_;
  ::xsd::cxx::tree::one< z_type > z_;
};

class initial_velocity: public ::xml_schema::type
{
  public:
  // x
  //
  typedef ::xml_schema::double_ x_type;
  typedef ::xsd::cxx::tree::traits< x_type, char, ::xsd::cxx::tree::schema_type::double_ > x_traits;

  const x_type&
  x () const;

  x_type&
  x ();

  void
  x (const x_type& x);

  // y
  //
  typedef ::xml_schema::double_ y_type;
  typedef ::xsd::cxx::tree::traits< y_type, char, ::xsd::cxx::tree::schema_type::double_ > y_traits;

  const y_type&
  y () const;

  y_type&
  y ();

  void
  y (const y_type& x);

  // z
  //
  typedef ::xml_schema::double_ z_type;
  typedef ::xsd::cxx::tree::traits< z_type, char, ::xsd::cxx::tree::schema_type::double_ > z_traits;

  const z_type&
  z () const;

  z_type&
  z ();

  void
  z (const z_type& x);

  // Constructors.
  //
  initial_velocity (const x_type&,
                    const y_type&,
                    const z_type&);

  initial_velocity (const ::xercesc::DOMElement& e,
                    ::xml_schema::flags f = 0,
                    ::xml_schema::container* c = 0);

  initial_velocity (const initial_velocity& x,
                    ::xml_schema::flags f = 0,
                    ::xml_schema::container* c = 0);

  virtual initial_velocity*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  initial_velocity&
  operator= (const initial_velocity& x);

  virtual 
  ~initial_velocity ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< x_type > x_;
  ::xsd::cxx::tree::one< y_type > y_;
  ::xsd::cxx::tree::one< z_type > z_;
};

class coordinate: public ::xml_schema::type
{
  public:
  // x
  //
  typedef ::xml_schema::double_ x_type;
  typedef ::xsd::cxx::tree::traits< x_type, char, ::xsd::cxx::tree::schema_type::double_ > x_traits;

  const x_type&
  x () const;

  x_type&
  x ();

  void
  x (const x_type& x);

  // y
  //
  typedef ::xml_schema::double_ y_type;
  typedef ::xsd::cxx::tree::traits< y_type, char, ::xsd::cxx::tree::schema_type::double_ > y_traits;

  const y_type&
  y () const;

  y_type&
  y ();

  void
  y (const y_type& x);

  // z
  //
  typedef ::xml_schema::double_ z_type;
  typedef ::xsd::cxx::tree::traits< z_type, char, ::xsd::cxx::tree::schema_type::double_ > z_traits;

  const z_type&
  z () const;

  z_type&
  z ();

  void
  z (const z_type& x);

  // Constructors.
  //
  coordinate (const x_type&,
              const y_type&,
              const z_type&);

  coordinate (const ::xercesc::DOMElement& e,
              ::xml_schema::flags f = 0,
              ::xml_schema::container* c = 0);

  coordinate (const coordinate& x,
              ::xml_schema::flags f = 0,
              ::xml_schema::container* c = 0);

  virtual coordinate*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  coordinate&
  operator= (const coordinate& x);

  virtual 
  ~coordinate ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< x_type > x_;
  ::xsd::cxx::tree::one< y_type > y_;
  ::xsd::cxx::tree::one< z_type > z_;
};

class dimensions: public ::xml_schema::type
{
  public:
  // x
  //
  typedef ::xml_schema::unsigned_int x_type;
  typedef ::xsd::cxx::tree::traits< x_type, char > x_traits;

  const x_type&
  x () const;

  x_type&
  x ();

  void
  x (const x_type& x);

  // y
  //
  typedef ::xml_schema::unsigned_int y_type;
  typedef ::xsd::cxx::tree::traits< y_type, char > y_traits;

  const y_type&
  y () const;

  y_type&
  y ();

  void
  y (const y_type& x);

  // z
  //
  typedef ::xml_schema::unsigned_int z_type;
  typedef ::xsd::cxx::tree::traits< z_type, char > z_traits;

  const z_type&
  z () const;

  z_type&
  z ();

  void
  z (const z_type& x);

  // Constructors.
  //
  dimensions (const x_type&,
              const y_type&,
              const z_type&);

  dimensions (const ::xercesc::DOMElement& e,
              ::xml_schema::flags f = 0,
              ::xml_schema::container* c = 0);

  dimensions (const dimensions& x,
              ::xml_schema::flags f = 0,
              ::xml_schema::container* c = 0);

  virtual dimensions*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  dimensions&
  operator= (const dimensions& x);

  virtual 
  ~dimensions ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< x_type > x_;
  ::xsd::cxx::tree::one< y_type > y_;
  ::xsd::cxx::tree::one< z_type > z_;
};

#include <iosfwd>

#include <xercesc/sax/InputSource.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMErrorHandler.hpp>

// Parse a URI or a local file.
//

::std::auto_ptr< ::MolSim >
MolSim_ (const ::std::string& uri,
         ::xml_schema::flags f = 0,
         const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::MolSim >
MolSim_ (const ::std::string& uri,
         ::xml_schema::error_handler& eh,
         ::xml_schema::flags f = 0,
         const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::MolSim >
MolSim_ (const ::std::string& uri,
         ::xercesc::DOMErrorHandler& eh,
         ::xml_schema::flags f = 0,
         const ::xml_schema::properties& p = ::xml_schema::properties ());

// Parse std::istream.
//

::std::auto_ptr< ::MolSim >
MolSim_ (::std::istream& is,
         ::xml_schema::flags f = 0,
         const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::MolSim >
MolSim_ (::std::istream& is,
         ::xml_schema::error_handler& eh,
         ::xml_schema::flags f = 0,
         const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::MolSim >
MolSim_ (::std::istream& is,
         ::xercesc::DOMErrorHandler& eh,
         ::xml_schema::flags f = 0,
         const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::MolSim >
MolSim_ (::std::istream& is,
         const ::std::string& id,
         ::xml_schema::flags f = 0,
         const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::MolSim >
MolSim_ (::std::istream& is,
         const ::std::string& id,
         ::xml_schema::error_handler& eh,
         ::xml_schema::flags f = 0,
         const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::MolSim >
MolSim_ (::std::istream& is,
         const ::std::string& id,
         ::xercesc::DOMErrorHandler& eh,
         ::xml_schema::flags f = 0,
         const ::xml_schema::properties& p = ::xml_schema::properties ());

// Parse xercesc::InputSource.
//

::std::auto_ptr< ::MolSim >
MolSim_ (::xercesc::InputSource& is,
         ::xml_schema::flags f = 0,
         const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::MolSim >
MolSim_ (::xercesc::InputSource& is,
         ::xml_schema::error_handler& eh,
         ::xml_schema::flags f = 0,
         const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::MolSim >
MolSim_ (::xercesc::InputSource& is,
         ::xercesc::DOMErrorHandler& eh,
         ::xml_schema::flags f = 0,
         const ::xml_schema::properties& p = ::xml_schema::properties ());

// Parse xercesc::DOMDocument.
//

::std::auto_ptr< ::MolSim >
MolSim_ (const ::xercesc::DOMDocument& d,
         ::xml_schema::flags f = 0,
         const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::MolSim >
MolSim_ (::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d,
         ::xml_schema::flags f = 0,
         const ::xml_schema::properties& p = ::xml_schema::properties ());

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.

#endif // MOL_SIM_HXX
