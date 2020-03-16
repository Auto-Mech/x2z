#include "chem.hh"
#include "units.hh"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <set>
#include <map>
#include <sstream>

/*********************** Atomic coordinates accuracies *******************/

double angle_tolerance = 5.0;

double distance_tolerance = 0.05;

std::set<std::set<int> > incipient_bond;

bool are_angles_equal (double a1, double a2)
{
  double da = a2 - a1;
  
  if(da < angle_tolerance && da > -angle_tolerance) {
    //
    return true;
  }
  
  return false;
}

bool are_distances_equal (double a1, double a2)
{
  double da = a2 - a1;
  
  if(da < distance_tolerance && da > -distance_tolerance) {
    //
    return true;
  }
  
  return false;
}

double max_bond_length(const AtomBase& a1, const AtomBase& a2)
{
  double res = atomic_radius(a1) + atomic_radius(a2);

  return 1.2 * res / Phys_const::bohr;
}

void MolecGeom::operator *= (const D3::Matrix& r)
{
  for(std::vector<Atom>::iterator at = begin(); at != end(); ++at)
    *at *= r;
}

void MolecGeom::operator += (const D3::Vector& d)
{
  for(std::vector<Atom>::iterator at = begin(); at != end(); ++at)
    *at += d;
}

void MolecGeom::operator -= (const D3::Vector& d)
{
  for(std::vector<Atom>::iterator at = begin(); at != end(); ++at)
    *at -= d;
}

void MolecGeom::operator *= (double d)
{
  for(std::vector<Atom>::iterator at = begin(); at != end(); ++at)
    *at *= d;
}

MolecOrient::MolecOrient (const MolecGeom& m)  : MolecGeom(m) 
{
  const char funame [] = "MolecOrient::MolecOrient (const MolecGeom&): ";

  int    itemp;
  
  double dtemp;
  
  bool   btemp;

  int min_ind;
  
  double min_val;
  bool istart;

  MolecGeom& m1 = *this;

  for(int i = 0; i < size(); ++i)
    if(m1[i] == AtomBase::DUMMY) {
      std::cout << funame << "no dummies, please\n";
      throw Error::Range();
    }

  if(size() < 2 ) {
    std::cout << funame << "wrong number of atoms\n";
    throw Error::Range();
  }

  m1 -= m[0];

  const double len0 = m1[1].vlength();
  if(are_distances_equal(len0, 0.)) {
    std::cout << funame << "atoms are too close\n";
    throw Error::Range();
  }

  if(size() == 2) {// diatomic
    _mt = LINEAR;
    for(int i = 0; i < 3; ++i)
      if(!i)
	m1[1][i] = len0;
      else
	m1[1][i] = 0.;
    return;
  }

  // is geometry linear?
  istart = true;
  for(int at = 2; at < size(); ++at) {
    dtemp = std::fabs(90. - angle(m1[1], m1[0], m1[at]));
    if(istart || dtemp < min_val) {
      istart = false;
      min_ind = at;
      min_val = dtemp;
    }
  }

  /************************ Linear geometry ******************************/

  if(are_angles_equal(min_val, 90.)) {// linear geometry
    _mt = LINEAR;

    // orient along x axis
    D3::Vector n0 = m1[1];
    n0 /= len0;
    for(int at = 1; at < size(); ++at)
      for(int i = 0; i < 3; ++i)
	if(!i)
	  m1[at][i] = vdot(m1[at], n0);
	else
	  m1[at][i] = 0.;

    // sorting atoms
    std::vector<double> coor(size());
    for(int at = 0; at < size(); ++at)
      coor[at] = m1[at][0];

    std::vector<int> perm(size());
    my_sort(coor, perm);

    MolecGeom m(size());
    
    for(int at = 0; at < size(); ++at)
      //
      m[at] = m1[perm[at]];
    
    m1 = m;

    /*
      std::cout << funame << "x-coordinates:";

    for(int a = 0; a < size(); ++a)
      //
      std::cout << "   " << m1[a][0];

    std::cout << "\n";
    */
    
    return;
  }
		
  /************************* Nonlinear geometry **************************/

  // change atoms order
  //
  if(min_ind != 2)
    //
    std::swap(m1[2], m1[min_ind]);

  m1 *= D3::Matrix(m1[1], m1[2]); // standard orientation

  // reference atoms distance matrix
  //
  for(int i = 0; i < 3; ++i)
    //
    _dm[i] = vdistance(m1[(i+1)%3], m1[(i+2)%3]);

  // is molecule plane?
  //
  _mt = PLANE;
  
  if(size() == 3)
    //
    return;

  for(int a = 3; a < size(); ++a)
    //
    if(!are_distances_equal(m1[a][2], 0.)) {
      //
      _mt = NONLINEAR;
      
      break;
    }
} 

int  MolecOrient::sym_num () const 
{
  return compare(*this, *this, SYMNUM); 
}

bool MolecOrient::is_enantiomer () const 
{
  if(_mt == NONLINEAR) {
    //
    MolecGeom m = *this;
    
    m *= -1.;
    
    return !compare(*this, MolecOrient(m), TEST);
  }
  
  return false;
}

int compare (const MolecOrient& m1, const MolecOrient& m2, 
	     MolecOrient::mode_t mode) 
{
  const char funame [] = "MolecOrient::compare: ";

  int    itemp;
  double dtemp;
  bool   btemp;

  if(m1.size() < 2 || m2.size() < 2) {
    //
    std::cerr << funame << "wrong number of atoms\n";
    
    throw Error::Range();
  }

  if(m1.size() != m2.size())
    //
    return 0;
    
  if(m1._mt != m2._mt)
    //
    return 0;

  int result = 0;

  /************************ Linear geometry ************************/

  if(m1._mt == MolecOrient::LINEAR) {
    //
    btemp = true;
    
    for(int at = 1; at < m1.size(); ++at)
      //
      if(m1[at] != m2[at] || !are_distances_equal(m1[at][0] - m1[0][0], m2[at][0] - m2[0][0])) {
	//
	btemp = false;

	break;
      }
    
    if(btemp) {
      //
      ++result;
      
      if(mode == MolecOrient::TEST)
	//
	return 1;
    }

    btemp = true;
    
    for(int at = 1; at < m1.size(); ++at)
      //
      if(m1[at] != m2[m1.size() - at - 1] || !are_distances_equal(m1[at][0] - m1[0][0], m2[m1.size() - 1][0] - m2[m1.size() - at - 1][0])) {
	//
	btemp = false;
	
	break;
      }
    
    if(btemp)
      //
      ++result;

    return result;
  }
    
  /************************ Nonlinear geometry ************************/

  //distance matrix
  //
  ConMat<double> tdm(m2.size());
  
  for(int i = 0; i < m2.size(); ++i)
    //
    for(int j = i+1; j < m2.size(); ++j)
      //
      tdm(i,j) = vdistance(m2[i], m2[j]);

  // three reference atoms cycle
  //
  int count = 0;
  
  for(int at0 = 0; at0 < m2.size(); ++at0) {
    //
    for(int at1 = 0; at1 < m2.size(); ++at1) {
      //
      if(at0 == at1)
	//
	continue;
      
      for(int at2 = 0; at2 < m2.size(); ++at2) {
	//
	if(at2 == at1 || at2 == at0)
	  //
	  continue;

	std::vector<int> perm(3);
	
	perm[0] = at0;
	perm[1] = at1;
	perm[2] = at2;

	// checking if reference atoms are the same
	//
	btemp = false;
	
	for(int i = 0; i < 3; ++i) {
	  //
	  if(m1[i] != m2[perm[i]]) {
	    //
	    btemp = true;
	    
	    break;
	  }
	}
	
	if(btemp)
	  //
	  continue;

	// checking if the distances between reference atoms are the same
	//
	double dm2[3];
	
	for(int i = 0; i < 3; ++i)
	  //
	  dm2[i] = tdm(perm[(i + 1) % 3], perm[(i + 2) % 3]);

	btemp = false;
	//
	for(int i = 0; i < 3; ++i) {
	  //
	  if(!are_distances_equal(m1._dm[i], dm2[i])) {
	    //
	    btemp = true;
	    
	    break;
	  }
	}
	
	if(btemp)
	  //
	  continue;

	// standard orientation
	//
	MolecGeom nm = m2;

	nm -= m2[perm[0]];

	nm *= D3::Matrix(nm[perm[1]], nm[perm[2]]);

	/*
	std::cout << "\n";
	
	for(int i = 0; i < 3; ++i)
	  //
	  std::cout << std::setw(3) << perm[i];
	
	std::cout << "\n";
	  
	for(int a = 0; a < m1.size(); ++a)
	  //
	  std::cout << m1[a] << "\n";
	
	std::cout << "\n";
		
	for(int a = 0; a < nm.size(); ++a)
	  //
	  std::cout << nm[a] << "\n";
	
	std::cout << "\n";
		
	*/
	
	// checking if the rest of the atoms coincide
	//
	for(int rest = 3; rest < m1.size(); ++rest) {
	  //
	  int    match = -1;
	  
	  double min_dist;
  
	  // best matching atom
	  //
	  for(int test = 0; test < nm.size(); ++test) {
	    //
	    btemp = false;
	    
	    for(int i = 0; i < perm.size(); ++i) {
	      //
	      if(test == perm[i]) {
		//
		btemp = true;
		
		break;
	      }
	    }
	    
	    if(btemp)
	      //
	      continue;
	    
	    if(m1[rest] != nm[test])
	      //
	      continue;
	    
	    dtemp = vdistance(m1[rest], nm[test]);
	    
	    if(match < 0 || dtemp < min_dist) {
	      //
	      min_dist = dtemp;
	      
	      match = test;
	    }
	    //
	  }// best matching atom
 
	  // different stoichiometry
	  //
	  if(match < 0) {
	    //
	    std::cerr << funame << "different stoichiometry\n";
	    
	    return 0;
	  }

	  // geometries differ
	  //
	  if(!are_distances_equal(min_dist, 0.)) {
	    //
	    //std::cout << "atom = " << rest << " best match = " << match << " distance = " << min_dist << "\n\n";
	    
	    break;
	  }
	  
	  perm.push_back(match);
	  //
	} // checking the rest of atoms

	// geometries are identical
	//
	if(perm.size() == m1.size()) {
	  //
	  if(mode == MolecOrient::TEST) {
	    //
	    return 1;
	  }
	  else
	    //
	    ++result;
	}
	
      }//
      //
    }//
    //
  }// reference atoms cycle
  
  return result;
}

PrimStruct::PrimStruct (const MolecGeom& g)
  //
  : ConMat<unsigned>(g.size()), MolecGeom(g), _la(g.size(), false)
{
  const char funame [] = "PrimStruct::PrimStruct(const MolecGeom&): ";

  typedef std::vector<Atom>::iterator mit_t;
  
  // check that there are no dummies
  //
  for(mit_t at = std::vector<Atom>::begin(); at != std::vector<Atom>::end(); ++at)
    //
    if(*at == AtomBase::DUMMY) {
      //
      std::cerr << funame << "WARNING: dummy atom encountered - removing\n";
      
      at = std::vector<Atom>::erase(at);
    }
  
  for(int i = 0; i < size(); ++i) {
    //
    for(int j = 0; j < i; ++j) {
      //
      // bond
      //
      std::set<int> bond;

      bond.insert(i);
      bond.insert(j);
      
      if(vdistance(g[i], g[j]) < max_bond_length(g[i], g[j]) || incipient_bond.find(bond) != incipient_bond.end()) {
	//
	(*this)(i, j) = 1;
      }
      else {
	//
	(*this)(i, j) = 0;
      }
    }
  }

  // check that the molecule is connected
  //
  //if(!is_connected())
    //
  //std::cerr << funame << "WARNING: bonding configuration is not connected\n";

  for(int i = 0; i < size(); ++i) {
    //
    // check that hydrogens are single bonded
    //
    if((*this)[i] == AtomBase::HYDROGEN && row_sum(i) > 1) {
      //
      std::cerr << funame << "WARNING: " << i << "-th hydrogen has more than one connection\n";
    }

    // check the number of bonds
    //
    if(row_sum(i) > 4) {
      //
      std::cerr << funame << "WARNING: " << i << "-th atom has more than four connections\n";
    }
    
    // check that oxygens have <= 2 connections
    //
    if((*this)[i] == AtomBase::OXYGEN && row_sum(i) > 2) {
      //
      std::cerr << funame << "WARNING: " << i << "-th oxygen has more than two connections\n";
    }
  }

  // check linearity
  //
  for(int at = 0; at < size(); ++at) {
    //
    // nearest neighbors
    
    std::vector<int> neighbor;
    //
    for(int at1 = 0; at1 < size(); ++at1) {
      //
      if(at1 != at && (*this)(at1, at))
	//
	neighbor.push_back(at1);
    }
    
    if(neighbor.size() == 2 && are_angles_equal(180., angle(g[neighbor[0]], g[at], g[neighbor[1]]))) {
      //
      _la[at] = true;

      //std::cout << funame << at << "-the atom is linear\n";
    }
  }
}

// check if the structure is connected
//
bool PrimStruct::is_connected (int at, const std::list<int>& group) const
{
  for(std::list<int>::const_iterator git = group.begin(); git != group.end(); ++git)
    //
    if(at == *git || is_connected(at, *git))
      //
      return true;

  return false;
}

std::list<std::list<int> > PrimStruct::connected_group () const
{
  int    itemp;
  bool   btemp;
  
  std::list<std::list<int> > res;
  
  std::list<int> pool;
  //
  for(int i = 0; i < size(); ++i)
    //
    pool.push_back(i);

  while(pool.size()) {
    //
    std::list<int> group;

    group.push_back(pool.back());

    pool.pop_back();
    
    btemp = true;

    while(btemp) {
      //
      btemp = false;

      for(std::list<int>::iterator pit = pool.begin(); pit != pool.end(); ++pit) {
	//
	if(is_connected(*pit, group)) {
	  //
	  btemp = true;
	
	  group.push_back(*pit);
	  
	  pool.erase(pit);

	  break;
	}
      }
    }

    res.push_back(group);
  }

  return res;
}

// check if the bond belong to a ring structure
//
bool PrimStruct::is_ring (int at1, int at2) const 
{
  const char funame [] = "PrimStruct::is_ring: ";

  if(at1 == at2) {
    //
    std::cerr << funame << "the same atom: " << at1 << "\n";

    throw Error::General();
  }

  if(!is_connected(at1, at2)) {
    //
    std::cerr << funame << "not connected: " << at1 << ", " << at2 << "\n";
			   
    throw Error::General();
  }

						       
  PrimStruct test = *this;
  
  test(at1, at2) = 0;

  return test.is_connected();
}

std::string PrimStruct::group_stoicheometry (const std::list<int>& group) const
{
  std::map<std::string, int> st;

  for(std::list<int>::const_iterator git = group.begin(); git != group.end(); ++git)
    //
    ++st[atom_name(*git)];

  std::ostringstream res;
  //
  for(std::map<std::string, int>::const_iterator sit = st.begin(); sit != st.end(); ++sit)
    //
    res << sit->first << sit->second;

  return res.str();
}
  
/*
// minimal number of bonds between two atoms
//
int  PrimStruct::distance (int at1, int at2) const 
{
  const char funame [] = "PrimStruct::distance: ";

  if(size() < 2)
    return 0;
  if(at1 == at2)
    return 0;

  std::vector<int> curr_set;
  curr_set.reserve(size());
  std::vector<int>  new_set;
  new_set.reserve(size());
  new_set.push_back(at1);

  int dist = 0;
  while(curr_set.size() != new_set.size()) {
    ++dist;
    curr_set = new_set;
    new_set.clear();
    for(int at = 0; at < size(); ++at)
      for(int i = 0; i < curr_set.size(); ++i)
	if(at == curr_set[i] || (*this)(at, curr_set[i])) {
	  if(at == at2) {
	    return dist;
	  }
	  new_set.push_back(at);
	  break;
	}
  }
  if(curr_set.size() != size()) {
    std::cout << funame << "primary structure is not connected\n";
    throw Error::Range();
  }
  std::cout << funame << "at2, " << at2 << ", is out of range\n";
  throw Error::Range();
}

*/

MolecStruct::MolecStruct (const PrimStruct& prim) 
   : PrimStruct(prim)
{
  const char funame [] = "MolecStruct::MolecStruct(const PrimStruct&): ";

  int    itemp;
  double dtemp;
  bool   btemp;
  
  ConMat<unsigned> m = prim;

  for(int i = 0; i < size(); ++i)
    //
    for(int j = 0; j < i; ++j)
      //
      if(m(i, j)) {
	//
	std::set<int> bond;

	bond.insert(i);
	
	bond.insert(j);

	if(incipient_bond.find(bond) == incipient_bond.end())
	  //
	  m(i, j) = 2;
      }

  // checking if the primary structure is sound
  //
  if(!is_connected()) {
    //
    std::cout << funame << "primary structure is not connected\n";
    
    throw Error::Range();
  }

  for(int i = 0; i < size(); ++i) {
    //
    if(m.row_sum(i) > 2 * valence(i)) {
      //
      std::cout << funame << i << "-th " << atom_name(i) << " atom real valence exceeds its formal valence\n";
      
      throw Error::Range();
    }
  }
  
  // getting all bonding configurations (resonances)
  //
  _resonance.push_back(m);

  // bond increment cycle
  //
  while(1) {
    //
    std::vector<ConMat<unsigned> > update;
    
    for(int r = 0; r < _resonance.size(); ++r)
      //
      for(int i = 0; i < size(); ++i)
	//
	for(int j = 0; j < i; ++j)
	  //
	  if((*this)(i, j) && _resonance[r].row_sum(i) + 1 < 2 * valence(i) && _resonance[r].row_sum(j) + 1 < 2 * valence(j)) {
	    //
	    update.push_back(_resonance[r]);
	    
	    (*update.rbegin())(i, j) += 2;
	  }
    
    // stop condition
    //
    if(!update.size())
      //
      break;

    // update resonance
    //
    _resonance.clear();
    //
    for(int r = 0; r < update.size(); ++r) {
      //
      bool btemp = true;
      
      for(int s = 0; s < _resonance.size(); ++s) {
	//
	if(_resonance[s] == update[r]) {
	  //
	  btemp = false;
	  
	  break;
	}
      }
      
      if(btemp)
	//
	_resonance.push_back(update[r]);
      //
    }//
    //
  }// bond increment cycle

  // first atom
  //
  _cpath.push_back(ConRec(0, -1));

  // atoms pool
  //
  std::vector<int> pool;
  
  for(int i = 1; i < size(); ++i)
    //
    pool.push_back(i);

  // reference (line #) of the previous connected atom
  //
  int cref = 0;

  // main cycle
  //
  while(pool.size()) {
    //
    // beginning of the new group
    //
    _cpath[cref].begin = _cpath.size();

    // previous atom
    //
    const int prev = _cpath[cref].atom;

    if(is_linear(prev)) {
      //
      _cpath[cref].attr |= LIN_BOND;
      
      _cpath.push_back(ConRec(-1, cref));
    }

    // pool cycle
    //
    for(std::vector<int>::iterator poolit = pool.begin(); poolit != pool.end(); ) {
      //
      // current atom
      //
      const int curr = *poolit;

      //new record
      //
      if(is_connected(curr, prev)) {
	//
	ConRec crec(curr, cref);
	
	// beta bond attribute
	//
	BetaData beta = is_beta(curr, prev);
	
	if(beta) {
	  //
	  _betvar[_cpath.size()] = beta;
	  
	  crec.attr |= BET_BOND;
	}

	_cpath.push_back(crec);
	
	poolit = pool.erase(poolit);
	//
      } // new record
      //
      else {
	//
	++poolit;
      }
      //
    } // pool cycle
    
    //
    // end of the new group
    //
    _cpath[cref].end = _cpath.size();

    while(_cpath[++cref].atom < 0) {}
    //
  } // main cycle

  // atom map
  for(int i = 0; i < _cpath.size(); ++i) {
    //
    itemp = _cpath[i].atom;
    
    if(itemp >= 0)
      //
      _atom_map[itemp] = i;
  }
    
  /*
  std::cout << "C_Path:\n"
	    << std::setw(7) << "line #"
	    << std::setw(7) << "atom"
	    << std::setw(7) << "cref"
	    << "\n";
  
  for(int i = 0; i < _cpath.size(); ++i)
    std::cout << std::setw(7) << i
	      << std::setw(7) << _cpath[i].atom
	      << std::setw(7) << _cpath[i].cref
	      << "\n";
  */
    
  _coval.resize(2, 3, _cpath.size());
  
  int  lroot = -1;

  bool lsingle;

  std::ostringstream to;
  //
  to << std::left;
  
  // constructing z-matrix
  //
  for(int ref0 = 0; ref0 < _cpath.size(); ++ref0) {
    //
    int ref1, ref2, ref3;
    
    // atom name
    //
    if(_cpath[ref0].atom < 0) {
      //
      to << "X ";
    }
    else {
      //
      to << std::setw(2) << atom_name(_cpath[ref0].atom);
    }
    
    // first reference (distance)
    //
    if(ref0 > 0) {
      //
      ref1 = _cpath[ref0].cref;
      
      to << ", " << std::setw(2) << ref1 + 1 << ", " << var_name(DISTANCE) << std::setw(2) << ref0;

      if(_cpath[ref0].atom < 0) {
	//
	_coval(DISTANCE, ref0) = 1.;

	_constvar.push_back(DISTANCE + 3 * ref0);
      }
      else {
	//
	_coval(DISTANCE, ref0) = ((*this)[_cpath[ref0].atom] - (*this)[_cpath[ref1].atom]).vlength();
      }
    }

    // second reference (polar angle)
    //
    if(ref0 > 1) {
      //
      // ref1 linear
      //
      if(_cpath[ref1].attr & LIN_BOND) {
	//
	// dummy atom
	//
	if(_cpath[ref0].atom < 0) {
	  //
	  ref2 = _cpath[ref1].cref;

	  _constvar.push_back(POLAR + 3 * ref0);
	}
	// real atom
	//
	else {
	  //
	  // dummy atom
	  //
	  ref2 = _cpath[ref1].begin;
	}
	
	_coval(POLAR, ref0) = 90.;
      }
      //
      // ref1 nonlinear
      //
      else {
	//
	// root atom
	//
	if(!ref1) {
	  //
	  ref2 = 1;
	}
	//
	else {
	  //
	  ref2 = _cpath[ref1].cref;
	}

	_coval(POLAR, ref0) = angle((*this)[_cpath[ref0].atom],
					  (*this)[_cpath[ref1].atom],
					  (*this)[_cpath[ref2].atom]);
      }
      
      to << ", " << std::setw(2) << ref2 + 1 << ", " << var_name(POLAR) << std::setw(2) << ref0;
    }

    // third reference (dihedral angle)
    //
    if(ref0 > 2) {
      //
      // ref1 linear
      //
      if(_cpath[ref1].attr & LIN_BOND) {
	//
	// dummy atom
	//
	if(_cpath[ref0].atom < 0) {
	  //
	  if(_cpath[ref2].attr & LIN_BOND) {
	    //
	    // ref2 dummy atom
	    //
	    ref3 = _cpath[ref2].begin;
	  }
	  else if(ref2) {
	    //
	    ref3 = _cpath[ref2].cref;
	  }
	  else if(ref1 == 1) {
	    //
	    ref3 = 2;
	  }
	  else {
	    //
 	    ref3 = 1;
	  }	  

	  _coval(DIHEDRAL, ref0) = 0.;

	  _constvar.push_back(DIHEDRAL + 3 * ref0);
	}
	//
	// real atom
	//
	else {
	  //
	  // ref1 is root
	  //
	  if(!ref1) {
	    //
	    ref3 = 2;
	  }
	  else {
	    //
	    ref3 = _cpath[ref1].cref;
	  }	 

	  _coval(DIHEDRAL, ref0) = 180.;
	}
      }
      //
      // ref1 nonlinear
      //
      else {
	//
	bool islin = false;

	bool isring = is_ring(_cpath[ref1].atom, _cpath[ref2].atom);
	
	bool isingle = is_single(_cpath[ref1].atom, _cpath[ref2].atom);

	//std::cout << "atoms: " << _cpath[ref1].atom << ", " << _cpath[ref2].atom << ": is ring? " << isring << std::endl;
	
	bool isrot = !isring && isingle;

	// not the first atom in the group
	//
	if(ref0 != _cpath[ref1].begin) {
	  //
	  if(!ref1) {
	    //
	    ref3 = 2;
	  }
	  else {
	    //
	    ref3 = _cpath[ref1].begin;
	  }

	  isrot = false;
	}
	// ref2 is linear
	//
	else if(_cpath[ref2].attr & LIN_BOND) {
	  //
	  islin = true;
	  
	  // ref2 dummy atom as dihedral reference ref3
	  //
	  ref3 = _cpath[ref2].begin;
	  
	  int curr = ref2;

	  int prev = ref1;

	  while(curr && _cpath[curr].attr & LIN_BOND) {
	    //
	    if(is_single(_cpath[curr].atom, _cpath[prev].atom))
	      //
	      isingle = true;

	    prev = curr;

	    curr = _cpath[curr].cref;
	  }
	  
	  isrot = !isring && isingle;

	  // nonlinear non-root atom
	  //
	  if(curr) {
	    //
	    curr = _cpath[curr].cref;

	    _coval(DIHEDRAL, ref0) = angle((*this)[_cpath[ref0].atom], 
						 (*this)[_cpath[ref1].atom],
						 (*this)[_cpath[ref2].atom], 
						 (*this)[_cpath[curr].atom]);
	  }
	  //
	  // root atom
	  //
	  else {
	    //
	    // linear root
	    //
	    if(_cpath[curr].attr & LIN_BOND) {
	      //
	      // first approach to the linear root
	      //
	      if(lroot < 0) {
		//
		lroot = _cpath[ref0].atom;

		lsingle = isingle;

		_coval(DIHEDRAL, ref0) = 0.;

		_constvar.push_back(DIHEDRAL + 3 * ref0);

		isrot = false;
	      }
	      // second approach to the linear root
	      //
	      else {
		//
		_coval(DIHEDRAL, ref0) = angle((*this)[_cpath[ref0].atom],
						     (*this)[_cpath[ref1].atom],
						     (*this)[_cpath[ref2].atom],
						     (*this)[lroot]);

		isrot = !isring && (isingle || lsingle);
	      }
	    }
	    // nonlinear root
	    //
	    else {
	      //
	      // terminal root
	      //
	      if(row_sum(0) == 1) {
		//
		_coval(DIHEDRAL, ref0) = 0.;

		_constvar.push_back(DIHEDRAL + 3 * ref0);

		isrot = false;
	      }
	      // root with several bonds
	      //
	      else if(prev == 1) {
		//
		_coval(DIHEDRAL, ref0) = angle((*this)[_cpath[ref0].atom], 
					       (*this)[_cpath[ref1].atom],
					       (*this)[_cpath[ref2].atom], 
					       (*this)[_cpath[2   ].atom]);
	      }
	      else {
		//
		_coval(DIHEDRAL, ref0) = angle((*this)[_cpath[ref0].atom],
					       (*this)[_cpath[ref1].atom],
					       (*this)[_cpath[ref2].atom],
					       (*this)[_cpath[1   ].atom]);
	      }
	    }
	  }
	  //
	  //
	} // ref2 is linear
	//
	else if(ref2) {
	  //
	  ref3 = _cpath[ref2].cref;
	}
	else if(ref1 == 1) {
	  //
	  ref3 = 2;
	}
	else {
	  //
	  ref3 = 1;
	}

	if(!islin) {
	  //
	  _coval(DIHEDRAL, ref0) = angle((*this)[_cpath[ref0].atom], 
					       (*this)[_cpath[ref1].atom],
					       (*this)[_cpath[ref2].atom], 
					       (*this)[_cpath[ref3].atom]);
	}

	if(isrot) {
	  //
	  PrimStruct test = *this;

	  test(_cpath[ref1].atom, _cpath[ref2].atom) = 0;
	  
	  _rotvar[ref0] = test.connected_group();

	  itemp = _rotvar[ref0].size();
	  
	  if(itemp != 2) {
	    //
	    std::cerr << funame << "wrong number of rotational groups: " << itemp << "\n";
	    
	    throw Error::General();
	  }
	  
	  for(std::list<std::list<int> >::iterator git = _rotvar[ref0].begin(); git != _rotvar[ref0].end(); ++git) {
	    //
	    for(std::list<int>::iterator it = git->begin(); it != git->end(); ++it)
	      //
	      if(*it == _cpath[ref1].atom || *it == _cpath[ref2].atom) {
		//
		itemp = *it;
		
		git->erase(it);
		
		git->push_front(itemp);
		
		break;
	      }
	  }
	}
	//
	//
      } // ref1 is nonlinear
      
      to << ", " << std::setw(2) << ref3 + 1 << ", " << var_name(DIHEDRAL) << std::setw(2) << ref0;
      //
      //
    } // dihedral angle value and reference

    to << "\n";
  }
  
  _zmat = to.str();
}

// atom-to-zmatrix map
//
int MolecStruct::atom_map (int i) const
{
  std::map<int, int>::const_iterator it = _atom_map.find(i);
  
  if(it != _atom_map.end())
    //
    return it->second;
  
  return -1;
}

void MolecStruct::print (std::ostream& to, const std::string& offset) const
{
  const PrimStruct& m = *this;
  
  double dtemp;
  bool   btemp;
  int    itemp;

  int old_precision = to.precision(2);

  to << offset << "Molecular structure:";

  if(_resonance.size() > 1) {
    //
    to << " resonantly stabilized (" << _resonance.size() << " resonances)" ;
  }

  to << "\n\n"; 

  to << offset << "   A\\A  ";

  for(int i = 0; i < size(); ++i) {
    //
    to << std::setw(3) <<  atom_name(i);

    if(i < 9) {
      //
      to << i + 1 << " ";
    }
    else
      //
      to << i + 1;
  }

  to << "\n\n";

  for(int j = 0; j < size(); ++j) {
    //
    to << offset << std::setw(4) << atom_name(j);

    if(j < 9) {
      //
      to << j + 1 << " ";
    }
    else
      //
      to << j + 1;

    for(int i = 0; i < j; ++i) {
      //
      to << std::setw(5) << bond_order(i, j);
    }

    to << std::setw(5) << "X" << "\n\n";
  }

  to << "\n";

  std::vector<int> rs_vec;

  for(int rad = 0; rad < size(); ++rad) {
    //
    btemp = false;

    for(int r = 0; r < _resonance.size(); ++r) {
      //
      if(_resonance[r].row_sum(rad) < 2 * valence(rad)) {
	//
	btemp = true;

	break;
      }
    }

    if(btemp)
      //
      rs_vec.push_back(rad);
  }

  if(rs_vec.size() > 1) {
    to << offset << "Free radical: Radical sites are ";
    for(int i = 0; i < rs_vec.size(); ++i) {
      itemp = rs_vec[i];
      to << atom_name(itemp) << itemp + 1 << " ";
    }
  }
  else if(rs_vec.size() == 1) {
    itemp = rs_vec[0];
    to << offset << "Free radical: Radical site is "
       << atom_name(itemp) << itemp + 1 << " ";
  }
  to << "\n\n";

  to << "Z-Matrix atom order:\n";

  for(int i = 0; i < _cpath.size(); ++i) {
    //
    to << std::setw(2) << i + 1 << " --> " << std::setw(2);

    if(_cpath[i].atom < 0) {
      //
      to << "X";
    }
    else
      //
      to << _cpath[i].atom + 1;
    
    to << "\n";
  }

  to << "\n";
  
  to.precision(old_precision);
}

std::vector<int> MolecStruct::atom_ordering() const
{
  unsigned int n = _cpath.size();
  std::vector<int> order(n);

  for(unsigned int i = 0; i < n; ++i)
    order[i] = _cpath[i].atom;

  return order;
}

// check if the bond is single
//
bool MolecStruct::is_single (int at0, int at1) const 
{
  const char funame [] = "MolecStruct::is_single: "; 

  if(!is_connected(at0, at1)) {
    //
    std::cerr << funame << "no bond";

    throw Error::Logic();
  }

  for(int r = 0; r < _resonance.size(); ++r)
    //
    if(_resonance[r](at0, at1) > 2)
      //
      return false;

  return true;
}

// check if the site is a radical one
bool MolecStruct::is_radical (int at) const
{
  for(int r = 0; r < _resonance.size(); ++r)
    //
    if(_resonance[r].row_sum(at) < 2 * valence(at))
      //
      return true;

  return false;
}

BetaData MolecStruct::is_beta (int at0, int at1) const 
{
  // is the bond single one
  //
  BetaData res;
  
  if(!is_single(at0, at1))
    //
    return res;

  // does the bond belong to a ring structure
  //
  if(is_ring(at0, at1))
    //
    res.isring = true;

  // find if there is a radical site next to the bond
  //
  for(int rad = 0; rad < size(); ++rad) {
    //
    if(rad == at0 || rad == at1)
      //
      continue;
    
    if(!is_radical(rad))
      //
      continue;
    
    if(is_connected(rad, at0)) {
      //
      res.radical = rad;

      res.primary = at0;

      res.secondary = at1;

      return res;
    }
    
    if(is_connected(rad, at1)) {
      //
      res.radical = rad;

      res.primary = at1;

      res.secondary = at0;

      return res;
    }
  }
  
  return res;
}

double MolecStruct::bond_order (int i, int j) const
{
  if(!(*this)(i,j))
    //
    return 0.;

  double res = 0.;

  for(int r = 0; r < _resonance.size(); ++r)
    //
    res += (double)_resonance[r](i, j);

  res /= (double)_resonance.size() * 2.;

  return res;
}

const char* MolecStruct::var_name (int v)
{
  const char funame [] = "MolecStruct::var_name: ";

  if(DISTANCE == v)
    return "R";

  if(POLAR == v)
    return "A";
  
  if(DIHEDRAL == v)
    return "D";
  
  std::cerr << funame << "out of range\n";
  throw Error::Range();
}

