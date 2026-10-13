// -*- C++ -*-

/*! @file
 * @brief Monomials - gauge action or fermion binlinear contributions for HMC
 */

#ifndef __abs_monomial_h__
#define __abs_monomial_h__

#include "wilstype_fermact_w.h"
#include "gaugeact.h"

#include "update/molecdyn/field_state.h"
#include "io/xmllog_io.h"
#include <cmath>

namespace Chroma
{
  //! Internal field refresh mode for pseudofermion updates
  enum InternalFieldsRefreshMode
  {
    INTERNAL_FIELDS_REFRESH_FULL,
    INTERNAL_FIELDS_REFRESH_OU
  };

  namespace MonomialInternalFieldsUtils
  {
    inline void appendDim(multi1d<int>& dims, int value)
    {
      multi1d<int> new_dims(dims.size() + 1);
      for (int i = 0; i < dims.size(); ++i) {
        new_dims[i] = dims[i];
      }
      new_dims[dims.size()] = value;
      dims = new_dims;
    }

    template<typename T>
    inline void captureDims(const T&, multi1d<int>&)
    {
      // Leaf type.
    }

    template<typename T>
    inline void captureDims(const multi1d<T>& fields, multi1d<int>& dims)
    {
      appendDim(dims, fields.size());
      if (fields.size() > 0) {
        captureDims(fields[0], dims);
      }
    }

    template<typename T>
    inline void resizeToDims(T&, const multi1d<int>& dims, int depth = 0)
    {
      if (depth != dims.size()) {
        throw std::string("MonomialInternalFieldsUtils::resizeToDims: bad depth for leaf type");
      }
    }

    template<typename T>
    inline void resizeToDims(multi1d<T>& fields, const multi1d<int>& dims, int depth = 0)
    {
      if (depth >= dims.size()) {
        throw std::string("MonomialInternalFieldsUtils::resizeToDims: missing dimension for multi1d type");
      }

      fields.resize(dims[depth]);
      for (int i = 0; i < fields.size(); ++i) {
        resizeToDims(fields[i], dims, depth + 1);
      }
    }

    template<typename T>
    inline void applyOU(T& refreshed_field, const T& old_field, const Real& c1, const Real& c2)
    {
      refreshed_field = c1 * old_field + c2 * refreshed_field;
    }

    template<typename T>
    inline void applyOU(multi1d<T>& refreshed_field, const multi1d<T>& old_field, const Real& c1, const Real& c2)
    {
      if (refreshed_field.size() != old_field.size()) {
        QDPIO::cerr << "MonomialInternalFieldsUtils::applyOU: incompatible field sizes" << std::endl;
        QDP_abort(1);
      }

      for (int i = 0; i < refreshed_field.size(); ++i) {
        applyOU(refreshed_field[i], old_field[i], c1, c2);
      }
    }

    inline void getOUCoeffs(const Real& traj_length, const Real& gamma, Real& c1, Real& c2)
    {
      c1 = exp(-gamma * traj_length);
      c2 = sqrt(Real(1) - c1 * c1);
    }

    template<typename T>
    inline void writeLeaves(QDPFileWriter& to, const T& field, int monomial_index, int& leaf_index)
    {
      XMLBufferWriter record_xml;
      push(record_xml, "PseudoFermionLeaf");
      write(record_xml, "MonomialIndex", monomial_index);
      write(record_xml, "LeafIndex", leaf_index);
      pop(record_xml);

      write(to, record_xml, field);
      ++leaf_index;
    }

    template<typename T>
    inline void writeLeaves(QDPFileWriter& to, const multi1d<T>& fields, int monomial_index, int& leaf_index)
    {
      for (int i = 0; i < fields.size(); ++i) {
        writeLeaves(to, fields[i], monomial_index, leaf_index);
      }
    }

    template<typename T>
    inline void readLeaves(QDPFileReader& from, T& field, int monomial_index, int& leaf_index)
    {
      XMLReader record_xml;
      read(from, record_xml, field);

      int file_monomial_index = -1;
      int file_leaf_index = -1;
      if (record_xml.count("/PseudoFermionLeaf/MonomialIndex") == 1) {
        read(record_xml, "/PseudoFermionLeaf/MonomialIndex", file_monomial_index);
      }
      if (record_xml.count("/PseudoFermionLeaf/LeafIndex") == 1) {
        read(record_xml, "/PseudoFermionLeaf/LeafIndex", file_leaf_index);
      }

      if (file_monomial_index != monomial_index || file_leaf_index != leaf_index) {
        throw std::string("MonomialInternalFieldsUtils::readLeaves: record index mismatch");
      }

      ++leaf_index;
    }

    template<typename T>
    inline void readLeaves(QDPFileReader& from, multi1d<T>& fields, int monomial_index, int& leaf_index)
    {
      for (int i = 0; i < fields.size(); ++i) {
        readLeaves(from, fields[i], monomial_index, leaf_index);
      }
    }

    template<typename T>
    inline void saveToOpenQIO(QDPFileWriter& to, const T& fields, int monomial_index)
    {
      multi1d<int> dims;
      captureDims(fields, dims);

      multi1d<int> meta(2 + dims.size());
      meta[0] = monomial_index;
      meta[1] = dims.size();
      for (int i = 0; i < dims.size(); ++i) {
        meta[2 + i] = dims[i];
      }

      XMLBufferWriter meta_xml;
      push(meta_xml, "PseudoFermionMeta");
      write(meta_xml, "MonomialIndex", monomial_index);
      pop(meta_xml);
      BinaryBufferWriter meta_bin;
      write(meta_bin, meta);
      write(to, meta_xml, meta_bin);

      int leaf_index = 0;
      writeLeaves(to, fields, monomial_index, leaf_index);
    }

    template<typename T>
    inline bool loadFromOpenQIO(QDPFileReader& from, T& fields, int expected_monomial_index)
    {
      XMLReader meta_xml;
      multi1d<int> meta;
      try {
        BinaryBufferReader meta_bin;
        read(from, meta_xml, meta_bin);
        read(meta_bin, meta);
      }
      catch (...) {
        return false;
      }

      if (meta.size() < 2) {
        return false;
      }

      const int monomial_index = meta[0];
      const int ndims = meta[1];
      if (monomial_index != expected_monomial_index || ndims < 0 || meta.size() != ndims + 2) {
        return false;
      }

      multi1d<int> dims(ndims);
      for (int i = 0; i < ndims; ++i) {
        dims[i] = meta[2 + i];
      }

      resizeToDims(fields, dims);

      int leaf_index = 0;
      try {
        readLeaves(from, fields, expected_monomial_index, leaf_index);
      }
      catch (...) {
        return false;
      }

      return true;
    }
  }

  //! An abstract monomial class, for inexact algorithms
  /*! @ingroup monomial
   *
   * Inexact in this case means energy computation is not supported,
   * (in an inexact algorithm sense -- obviously it is weird to have
   * a hamiltonian where you cannot compute the energy. We may need
   * to think more about this)
   *
   * This serves the following purpose. It definees
   * an interface for computing the total force
   * and can refresh the momenta,
   *
   *
   * We don't specify how the momenta is refreshed. It is "virtual".
   * HMD type algorithms will porbably use gaussian noise.
   * GHMD type algorithms will mix the previous momenta some. How
   * to do that will be encoded in the derived class, probably
   * through the constructor.
   *
   *
   * For this it needs to know the types of coordinates and the momenta
   * so that it can act on the right kind of state.
   */
  template<typename P, typename Q>
  class Monomial
  {
  public:
    //! virtual destructor:
    virtual ~Monomial() {}

    //! Compute dsdq for the system...
    /*! Not specified how to actually do this s is the state, F is the computed force */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermion fields if any
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Refresh pseudofermion fields with optional OU parameters
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state,
                                       const Real& traj_length,
                                       const Real& gamma,
                                       const InternalFieldsRefreshMode mode)
    {
      refreshInternalFields(field_state);
    }

    //! Copy pseudofermion fields from another monomial...
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

    //! Whether this monomial carries internal fields
    virtual bool hasInternalFields(void) const
    {
      return false;
    }

    //! Push current internal fields to a backup stack
    virtual void pushInternalFields(void) {}

    //! Pop internal fields from a backup stack (restore)
    virtual void popInternalFields(void) {}

    //! Drop top internal-field backup without restore
    virtual void dropInternalFields(void) {}

    //! Save internal fields to an open QIO stream
    virtual bool saveInternalFields(QDPFileWriter& to, int monomial_index) const
    {
      return false;
    }

    //! Load internal fields from an open QIO stream
    virtual bool loadInternalFields(QDPFileReader& from, int monomial_index)
    {
      return false;
    }

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }
  };


  //-------------------------------------------------------------------------------------------
  //! Abstract monomial class, for exact algorithms
  /*! @ingroup monomial
   *
   * Now define similar classes for exact algorithms.
   * These are basically the same as before but they can compute
   * energies too. Do these need to inerit?
   * Yes. Reason: We can always give it to an inexact algorithm through
   * a downcast. In that case the energy calculations would be hidden.
   */
  template<typename P, typename Q>
  class ExactMonomial : public Monomial<P, Q>
  {
  public:
    //! virtual destructor:
    virtual ~ExactMonomial() {}

    //! Compute dsdq for the system... Not specified how to actually do this
    //  s is the state, F is the computed force
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    // Compute the energies

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermion fields if any
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Refresh pseudofermion fields with optional OU parameters
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state,
                                       const Real& traj_length,
                                       const Real& gamma,
                                       const InternalFieldsRefreshMode mode)
    {
      Monomial<P,Q>::refreshInternalFields(field_state, traj_length, gamma, mode);
    }

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }
  };

  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup monomial
   *
   * The fermion hierachy would splits at the very top into
   * inexact and exact monomials. An exact monomial can be used
   * for an inexact algorithm, but not vice-versa.
   */

  /* Unfortunately we need to template on the Phi-s because
     we need that template for the FermActs */
  template<typename P, typename Q, typename Phi>
  class FermMonomial : public Monomial<P,Q>
  {
  public:
    //! virtual destructor:
    ~FermMonomial() {}

    //! Compute dsdq for the system... Not specified how to actually do this
    //  s is the state, F is the computed force
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) = 0;

    // Refresh all pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0 ;

    //! Refresh pseudofermion fields with optional OU parameters
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state,
                                       const Real& traj_length,
                                       const Real& gamma,
                                       const InternalFieldsRefreshMode mode)
    {
      Monomial<P,Q>::refreshInternalFields(field_state, traj_length, gamma, mode);
    }

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }
  };


  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup monomial
   *
   * The fermion hierachy would splits at the very top into
   * inexact and exact monomials. An exact monomial can be used
   * for an inexact algorithm, but not vice-versa.
   *
   * Unfortunately we need to template on the Phi-s because
   *  we need that template for the FermActs
   */
  template<typename P, typename Q, typename Phi>
  class ExactFermMonomial : public ExactMonomial<P,Q>
  {
  public:
    //! virtual destructor:
    ~ExactFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    /*! s is the state, F is the computed force */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Refresh pseudofermion fields with optional OU parameters
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state,
                                       const Real& traj_length,
                                       const Real& gamma,
                                       const InternalFieldsRefreshMode mode)
    {
      Monomial<P,Q>::refreshInternalFields(field_state, traj_length, gamma, mode);
    }

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }
  };


  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup monomial
   *
   * Exact fermionic monomials with pseudofermions living in 4D
   *
   * We need to template on the Phi-s because of the fermacts
   */
  template<typename P, typename Q, typename Phi>
  class ExactFermMonomial4D : public ExactFermMonomial<P,Q,Phi>
  {
  public:
    //! virtual destructor:
    ~ExactFermMonomial4D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    /*! s is the state, F is the computed force */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Refresh pseudofermion fields with optional OU parameters
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state,
                                       const Real& traj_length,
                                       const Real& gamma,
                                       const InternalFieldsRefreshMode mode)
    {
      Monomial<P,Q>::refreshInternalFields(field_state, traj_length, gamma, mode);
    }

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }
  };


  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup monomial
   *
   * Exact fermionic monomials with pseudofermions living in 4D
   *
   * We need to template on the Phi-s because of the fermacts
   */
  template<typename P, typename Q, typename Phi>
  class ExactFermMonomial5D : public ExactFermMonomial<P,Q,Phi>
  {
  public:
    //! virtual destructor:
    ~ExactFermMonomial5D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    //  s is the state, F is the computed force
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Refresh pseudofermion fields with optional OU parameters
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state,
                                       const Real& traj_length,
                                       const Real& gamma,
                                       const InternalFieldsRefreshMode mode)
    {
      Monomial<P,Q>::refreshInternalFields(field_state, traj_length, gamma, mode);
    }

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }
  };


  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup monomial
   *
   * Wilson-like fermion monomials. Not sure what these really do that
   * is new. There can be a staggered version.
   */
  template<typename P, typename Q, typename Phi>
  class ExactWilsonTypeFermMonomial : public ExactFermMonomial4D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~ExactWilsonTypeFermMonomial() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    /*! s is the state, F is the computed force */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Refresh pseudofermion fields with optional OU parameters
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state,
                                       const Real& traj_length,
                                       const Real& gamma,
                                       const InternalFieldsRefreshMode mode)
    {
      Monomial<P,Q>::refreshInternalFields(field_state, traj_length, gamma, mode);
    }

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }

  protected:
    //! Get at fermion action for pseudofermion field i
    virtual const WilsonTypeFermAct<Phi,P,Q>& getFermAct(void) const = 0;

  };


  //-------------------------------------------------------------------------------------------
  //! Fermionic monomials (binlinears in fermion fields)
  /*! @ingroup monomial
   *
   * Wilson-like fermion monomials. Not sure what these really do that
   * is new. There can be a staggered version.
   */
  template<typename P, typename Q, typename Phi>
  class ExactWilsonTypeFermMonomial5D : public ExactFermMonomial5D<P,Q,Phi>
  {
  public:
     //! virtual destructor:
    ~ExactWilsonTypeFermMonomial5D() {}

    //! Compute the total action
    virtual Double S(const AbsFieldState<P,Q>& s)  = 0;

    //! Compute dsdq for the system... Not specified how to actually do this
    /*! s is the state, F is the computed force */
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s)  = 0;

    //! Refresh pseudofermions
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;

    //! Refresh pseudofermion fields with optional OU parameters
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state,
                                       const Real& traj_length,
                                       const Real& gamma,
                                       const InternalFieldsRefreshMode mode)
    {
      Monomial<P,Q>::refreshInternalFields(field_state, traj_length, gamma, mode);
    }

    //! Copy pseudofermions if any
    virtual void setInternalFields(const Monomial<P,Q>& m) = 0;

    //! Reset predictors
    virtual void resetPredictors(void) { /* Nop for most */ }

  protected:
    //! Get at fermion action for pseudofermion field i
    virtual const WilsonTypeFermAct5D<Phi,P,Q>& getFermAct(void) const = 0;

  };

}


#endif
