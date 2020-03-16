#include<vector>
#include <string>
#include<fstream>
#include<iostream>
#include<sstream>
#include<set>

#include "libx2z/units.hh"
#include "libx2z/chem.hh"
#include "libx2z/math.hh"

int main(int argc, const char* argv [])
{
  const char funame [] = "x2z: ";

  int    itemp;
  double dtemp;
  
  if(argc != 2) {
    //
    std::cout << "usage: x2z input_file\n";
    
    return 1;
  }
  
  std::ifstream from(argv[1]);
  
  if(!from) {
    //
    std::cout << funame << "cannot open " << argv[1] << " file\n";

    return 1;
  }

  std::string token, line, comment;

  // read molecular geometry
  //
  MolecGeom geom;

  if(!(from >> itemp)) {
    //
    std::cout << funame << "cannot read number of atoms\n";
    
    return 1;
  }
  
  std::getline(from, comment);

  std::getline(from, comment);

  for(int i = 0; i < itemp; ++i) {
    //
    Atom a(from);

    a /= Phys_const::bohr;
    
    geom.push_back(a);
  }
  
  if(!from) {
    //
    std::cout << funame << " cannot read molecular geometry\n";

    return 1;
  }

  // auxiliary input
  //
  std::set<std::string> avail_keys;
  
  const std::string atol_key = "AngleTolerance";
  const std::string bond_key = "IncipientBond";

  avail_keys.insert(atol_key);
  avail_keys.insert(bond_key);

  // input cycle
  //
  while(from >> token) {
    //
    // angle tolerance
    //
    if(atol_key == token) {
      //
      if(!(from >> dtemp)) {
	//
	std::cerr << funame << token << ": corrupted\n";
	
	return 1;
      }
      
      if(dtemp <= 0.) {
	//
	std::cerr << funame << token << ": out of range\n";

	return 1;
      }

      angle_tolerance = dtemp;
    }
    // incipient bond
    //
    else if(bond_key == token) {
      //
      std::set<int> bond;

      std::getline(from, line);

      std::istringstream iss(line);
      
      while(iss >> itemp) {
	//
	if(itemp < 1 || itemp > geom.size()) {
	  //
	  std::cerr << funame << token << ": out of range: " << itemp << "\n";

	  return 1;
	}

	if(!bond.insert(--itemp).second) {
	  //
	  std::cerr << funame << token << ": duplicated index: " << itemp + 1 << "\n";

	  return 1;
	}
      }

      if(bond.size() != 2) {
	//
	std::cerr << funame << token << ": wrong number of atoms: " << bond.size() << "\n";

	return 1;
      }

      incipient_bond.insert(bond);
    }
    // unknown keyword
    //
    else {
      //
      std::cerr << funame << "unknown keyword: " << token << "\nAvailable keywords:";


      for(std::set<std::string>::const_iterator kit = avail_keys.begin(); kit != avail_keys.end(); ++kit)
	//
	std::cerr << "   " << *kit;

      std::cerr << "\n";
      
      return 1;
    }
  }
  
  MolecOrient mo(geom);

  std::cout << "molecule is ";
  
  if(mo.is_linear()) {
    //
    std::cout << "linear\n";
  }
  else if(mo.is_plane()) {
    //
    std::cout << "plane\n";
  }
  else {
    //
    std::cout << "nonlinear\n";
    
    std::cout << "has enantiomer? ";
    
    if(mo.is_enantiomer()) {
      //
      std::cout << "yes\n";
    }
    else
      //
      std::cout << "no\n";
  }
  
  std::cout << "\n";

  std::cout << "rotational symmetry number = " << mo.sym_num() << "\n\n";

  PrimStruct prim(geom);

  if(!prim.is_connected()) {
    std::cout << funame << "primary structure is not connected\n";
    return 0;
  }

  MolecStruct mol(prim);

  mol.print(std::cout);

  std::cout << "Z-Matrix:\n";
  //
  std::cout << mol.zmatrix();
  //
  std::cout << "\n";

  std::cout << std::left;

  for(int i = 1; i < mol.zmat_coval().size(1); ++i) {
    //
    std::cout << MolecStruct::var_name(MolecStruct::DISTANCE) << std::setw(2) << i
	      << " = " << std::setw(9) << mol.zmat_coval()(MolecStruct::DISTANCE, i);
    
    if(i > 1) {
      //
      std::cout << " " << MolecStruct::var_name(MolecStruct::POLAR) <<  std::setw(2) << i
		<< " = " << std::setw(9) << mol.zmat_coval()(MolecStruct::POLAR, i);
    }

    if(i > 2) {
      //
      std::cout << " " <<  MolecStruct::var_name(MolecStruct::DIHEDRAL) << std::setw(2) << i
		<< " = " << std::setw(15) << mol.zmat_coval()(MolecStruct::DIHEDRAL, i);
    }

    std::cout << "\n";
  }
  std::cout << "\n";

  // constant parameters
  //
  if(mol.const_var().size()) {
    //
    std::cout << "Constants:";

    for(std::list<int>::const_iterator cit = mol.const_var().begin(); cit != mol.const_var().end(); ++cit) {
      //
      std::cout << "   " << MolecStruct::var_name(*cit % 3) << *cit / 3;
    }

    std::cout << "\n\n";
  }
  else
    std::cout << "no constants\n";

  // rotational dihedral angles
  //
  if(mol.rotation_bond().size()) {
    //
    std::cout << "Rotational bond dihedral angles: ";

    for(std::map<int, std::list<std::list<int> > >::const_iterator bit = mol.rotation_bond().begin();
	bit != mol.rotation_bond().end(); ++bit) {
      //
      if(bit != mol.rotation_bond().begin())
	//
	std::cout << ", ";

      std::cout << MolecStruct::var_name(MolecStruct::DIHEDRAL) << bit->first;
    }

    std::cout << "\n\n";

    std::cout << "Rotational groups:\n";

    for(std::map<int, std::list<std::list<int> > >::const_iterator bit = mol.rotation_bond().begin();
	bit != mol.rotation_bond().end(); ++bit) {
      
      std::cout << MolecStruct::var_name(MolecStruct::DIHEDRAL) << std::setw(2) << bit->first;

      for(std::list<std::list<int> >::const_iterator git = bit->second.begin(); git != bit->second.end(); ++git)
	//
	std::cout << " " << std::setw(10) << mol.group_stoicheometry(*git);

      std::cout << "\n";
    }

    std::cout << "\n";

    for(std::map<int, std::list<std::list<int> > >::const_iterator bit = mol.rotation_bond().begin();
	bit != mol.rotation_bond().end(); ++bit) {
      
      std::cout << MolecStruct::var_name(MolecStruct::DIHEDRAL) << std::setw(2) << bit->first;

      for(std::list<std::list<int> >::const_iterator git = bit->second.begin(); git != bit->second.end(); ++git) {
	//
	std::cout << "   ";
	
	for(std::list<int>::const_iterator it = git->begin(); it != git->end(); ++it) {
	  //
	  if(it != git->begin())
	    //
	    std::cout << ",";
	  
	  std::cout << mol[*it].name() << mol.atom_map(*it) + 1;
	}
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
  else
    std::cout << "no rotational bonds\n";
  
  // beta bonds
  //
  itemp = 0;
    
  for(std::map<int, BetaData>::const_iterator bit = mol.beta_bond().begin(); bit != mol.beta_bond().end(); ++bit) {
    //
    if(bit->second.isring)
      //
      continue;

    if(itemp++) {
      //
      std::cout << ", ";
    }
    else
      //
      std::cout << "Beta-scission bonds:   ";

      
    std::cout << MolecStruct::var_name(MolecStruct::DISTANCE) << bit->first;
  }

  if(itemp)
    //
    std::cout << "\n";

  
  return 0;

}
