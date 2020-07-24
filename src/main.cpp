#include <iostream>
#include "casm/global/definitions.hh"
#include <boost/filesystem.hpp>
#include "casm/crystallography/io/VaspIO.hh"

/// What is being tested:
#include "casm/crystallography/StrucMapping.hh"

/// What is being used to test it:
#include "casm/crystallography/Adapter.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/SimpleStrucMapCalculator.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/clex/io/json/ConfigMapping.hh"
#include "casm/clex/ConfigMapping.hh"

using namespace CASM;
using namespace xtal;

void print_mapping_node(xtal::MappingNode const &el, std::string label, std::ostream & stream = std::cout) {
  stream << "ELEMENT " << label << " ->  ";
  stream << "   cost: " << el.cost << "  bcost: " << el.atomic_node.cost << "  lcost: " << el.lattice_node.cost << "\n"
            << "   translation: " << el.atomic_node.translation.transpose() << "\n"
            << "   isometry: \n" << el.lattice_node.isometry << "\n"
            << "   stretch: \n" << el.lattice_node.stretch << "\n"
            << "   parent: \n" << el.lattice_node.parent.superlattice().lat_column_mat() << "\n"
            << "   child: \n" << el.lattice_node.child.superlattice().lat_column_mat() << "\n"
            << "   cost_mat: \n" << el.atomic_node.cost_mat << "\n"
            << "   partitioned: " << el.is_partitioned << "\n"
            << "   forced_on: \n";
  for(auto const &pr : el.atomic_node.forced_on)
    stream << "     (" << pr.first << ", " << pr.second << ")\n";
  stream << "   irow: " << el.atomic_node.irow << "\n"
            << "   icol: " << el.atomic_node.icol << "\n"
            << "   assignment: " << el.atomic_node.assignment << "\n"
            << "   displacement: \n" << el.atom_displacement << "\n"
            << "   tot assignment: " << el.atom_permutation << "\n"
            << "   MolMap: ";
  for(auto const &m : el.mol_map){
    stream << "{";
      for(auto const & a : m){
        stream << a << ".";
      }
      stream << "} ";
  }
  stream << "\n\n-----\n\n";

}

int main(){

  /**
   * This is where you get to shine. Place all your code in here
   * and use the infrastructure of the repository to compile
   * and install.
   *
   * If you're interested in working with something more complicated
   * than a single main.cpp file is convenient for, check out
   * casm-utilities!
   * https://github.com/goirijo/casm-utilities
   */
  
  ConfigMapping::Settings set;
  
  if(!fs::is_regular_file("mapping.json")){
    std::cerr<< "FILE REQUIRED: 'mapping.json'\n";
    return 1;
  }
  jsonParser kwargs;
  kwargs.read("mapping.json");

  std::string parent_path=kwargs["parent_path"].get<std::string>();
  xtal::BasicStructure parent;
  {
    std::ifstream pstream(parent_path);
    parent=BasicStructure::from_poscar_stream(pstream);
  }

  auto parent_fg = make_factor_group(parent);


  std::string child_path=kwargs["child_path"].get<std::string>();
  xtal::BasicStructure child;
  {
    std::ifstream cstream(child_path);
    child=BasicStructure::from_poscar_stream(cstream);
  }

  auto child_fg = make_factor_group(child);


  from_json(set,kwargs);
  
  //Number of mappings to find for each superlattice
  Index k_best=1;
  if(kwargs.contains("k_best")){
    k_best=kwargs["k_best"].get<Index>();
  }
  else{
    kwargs["k_best"]=k_best;
  }

  // Integer range of unimodular matrix
  Index trans_range=1;
  if(kwargs.contains("trans_range")){
    trans_range=kwargs["trans_range"].get<Index>();
  }
  else{
    kwargs["trans_range"]=trans_range;
  }
  
  std::string output_directory=parent_path+"."+child_path;
  for(char & c : output_directory){
    if(c=='/')
      c='.';
  }
  
  if(kwargs.contains("output_directory")){
    output_directory=kwargs["output_directory"].get<std::string>();
  }
  else{
    kwargs["output_directory"]="";
  }
  
  Index max_vol=1;
  if(kwargs.contains("max_vol")){
    max_vol=kwargs["max_vol"].get<Index>();
  }
  else{
    kwargs["max_vol"]=max_vol;
  }
  

  Index min_vol=1;
  if(kwargs.contains("min_vol")){
    min_vol=kwargs["min_vol"].get<Index>();
  }
  else{
    kwargs["min_vol"]=min_vol;
  }
  
  Index interp_images=0;
  if(kwargs.contains("interp_images")){
    interp_images=kwargs["interp_images"].get<Index>();
  }
  else{
    kwargs["interp_images"]=interp_images;
  }
  

  
  double min_cost(-1);
  double max_cost(1e10);
  if(kwargs.contains("max_cost")){
    max_cost=kwargs["max_cost"].get<double>();
  }
  else{
    kwargs["max_cost"]=max_cost;
  }
  
  if(kwargs.contains("min_cost")){
    min_cost=kwargs["min_cost"].get<double>();
  }
  else{
    kwargs["min_cost"]=min_cost;
  }

  bool primitive_only=true;
  if(kwargs.contains("primitive_only")){
    primitive_only=kwargs["primitive_only"].get<bool>();
  }
  else{
    kwargs["primitive_only"]=primitive_only;
  }

  bool write_files=true;
  if(kwargs.contains("write_files")){
    write_files=kwargs["write_files"].get<bool>();
  }
  else{
    kwargs["write_files"]=write_files;
  }

  bool use_parent_sym=true;
  if(kwargs.contains("use_parent_sym")){
    use_parent_sym=kwargs["use_parent_sym"].get<bool>();
  }
  else{
    kwargs["use_parent_sym"]=use_parent_sym;
  }

  bool use_child_sym=true;
  if(kwargs.contains("use_child_sym")){
    use_child_sym=kwargs["use_child_sym"].get<bool>();
  }
  else{
    kwargs["use_child_sym"]=use_child_sym;
  }
  

  
  bool midpoint_analysis=false;
  if(kwargs.contains("midpoint_analysis")){
    midpoint_analysis=kwargs["midpoint_analysis"].get<bool>();
  }
  else{
    kwargs["midpoint_analysis"]=midpoint_analysis;
  }
  

  std::vector<double> ca_range({1.,1e-6,1.});
  if(kwargs.contains("ca_range")){
    ca_range.clear();
    from_json(ca_range,kwargs["ca_range"]);
  }

  ca_range[1]=max(ca_range[1],1e-6);
  if(almost_equal(ca_range[0],ca_range[2],5e-7))
    ca_range[2]+=5e-7;
  
  kwargs["ca_range"]=ca_range;

  
  std::cout << "Default settings: \n";
  {
    jsonParser tjson;
    tjson=ConfigMapping::Settings();
    std::cout<< tjson << std::endl << std::endl;
  }

  std::cout << "Actual settings used: \n";
  std::cout<< kwargs << std::endl << std::endl;

  
  xtal::SimpleStructure pstruc = make_simple_structure(parent);

  xtal::SimpleStructure cstruc = make_simple_structure(child);

  int options=0;
  if(set.robust){
    options=StrucMapper::robust;
  }
  while(ca_range[0]<ca_range[2]){
    std::cout << "CA= " <<ca_range[0] << std::endl;
    xtal::StrucMapper mapper(xtal::SimpleStrucMapCalculator(pstruc,
                                                            use_parent_sym ? parent_fg : xtal::SymOpVector({xtal::SymOp::identity()}),
                                                            xtal::SimpleStructure::SpeciesMode::ATOM,
                                                            allowed_molecule_names(parent)),
                             set.lattice_weight,
                             set.max_vol_change,
                             options,
                             set.cost_tol,
                             set.min_va_frac,
                             set.max_va_frac);

    mapper.set_lattice_transformation_range(trans_range);

    //std::cout << "Number of maps: " << result.size() << std::endl;
    SuperlatticeEnumerator slat_enum(child_fg.begin(), child_fg.end(), child.lattice(), ScelEnumProps(min_vol, max_vol+1));
    Index i=0;
    Index s=0;
    
    cstruc.lat_column_mat = child.lattice().lat_column_mat()/sqrt(ca_range[0]);
    cstruc.lat_column_mat.col(2)=ca_range[0]*child.lattice().lat_column_mat().col(2);
    double real_ca=cstruc.lat_column_mat.col(2).norm()/cstruc.lat_column_mat.col(0).norm();
    std::cout << "Strained lattice volume: " << double(cstruc.lat_column_mat.determinant()) << std::endl;
    for(auto it=slat_enum.begin(); it!=slat_enum.end(); ++it){
      std::cout << "New slat\n";
      SimpleStructure super_cstruc=make_superstructure(it.matrix(),cstruc);

      auto tresult = mapper.map_deformed_struc(super_cstruc,
                                               k_best,
                                               max_cost,
                                               min_cost,
                                               false,
                                               use_child_sym?child_fg : decltype(child_fg){child_fg[0]});

      ++s;
    
      for(auto it2=tresult.begin(); it2!=tresult.end(); ++it2,++i){
        SimpleStructure mapped_cstruc=mapper.calculator().resolve_setting(*it2,super_cstruc);
        mapped_cstruc.atom_info=mapped_cstruc.mol_info;
        SimpleStructure super_pstruc=make_superstructure((it2->lattice_node).parent.transformation_matrix_to_super().cast<int>(),pstruc);
        SimpleStructure tstruc=super_pstruc;

        tstruc.lat_column_mat=0.1*super_pstruc.lat_column_mat+0.9*mapped_cstruc.lat_column_mat;
        tstruc.atom_info.coords=0.1*super_pstruc.atom_info.coords+0.9*mapped_cstruc.atom_info.coords;
        tstruc.mol_info.coords=0.1*super_pstruc.mol_info.coords+0.9*mapped_cstruc.mol_info.coords;

        bool is_invalid=primitive_only;
        if(primitive_only){
          xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(tstruc)));
          xtal::LatticeNode tnode((xtal::Lattice(tstruc.lat_column_mat)),
                                  (xtal::Lattice(tstruc.lat_column_mat)),
                                  (xtal::Lattice(tstruc.lat_column_mat)),
                                  (xtal::Lattice(tstruc.lat_column_mat)),
                                  tstruc.atom_info.size());
          auto trans_set = mapper.map_deformed_struc_impose_lattice_node(tstruc, tnode, 0, xtal::StrucMapping::big_inf(), 1e-3);
          if(trans_set.size()<=1){
            is_invalid=false;
          }
        }
        if(is_invalid)
          continue;
        
        std::stringstream ss;
        
        ss << i <<  "  vol: " << it.volume() << "  slat: " << s << "  " << "  CA: " << real_ca << "  ";
        print_mapping_node(*it2, ss.str());
        
        std::string map_dir=output_directory + "/" + std::to_string(i);
        if(write_files||midpoint_analysis || interp_images){
          fs::create_directories(map_dir);
          if(midpoint_analysis){
            tstruc.lat_column_mat=0.5*super_pstruc.lat_column_mat+0.5*mapped_cstruc.lat_column_mat;
            tstruc.atom_info.coords=0.5*super_pstruc.atom_info.coords+0.5*mapped_cstruc.atom_info.coords;
            tstruc.mol_info.coords=0.5*super_pstruc.mol_info.coords+0.5*mapped_cstruc.mol_info.coords;
            
            
            std::string struc_path=map_dir + "/midpoint.vasp";
            {
              std::ofstream fout(struc_path);
              VaspIO::PrintPOSCAR POS(tstruc);
              POS.print(fout);
            }
            
            xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(tstruc)));
            auto symresult = mapper.map_deformed_struc(tstruc,
                                                       0,
                                                       1.,
                                                       1.e-3,
                                                       false);
            {
              std::ofstream fout(map_dir + "/factor_group");
              for(auto const & el : symresult){
                print_mapping_node(el,"",fout);
              }
            }
          }
          
          

          SimpleStructure tstruc2=super_pstruc;
          tstruc2.atom_info.names.resize(2*tstruc.atom_info.names.size(),"H");
          tstruc2.mol_info.names.resize(2*tstruc.mol_info.names.size(),"H");
          tstruc2.atom_info.coords.resize(3,2*tstruc.atom_info.names.size());
          tstruc2.mol_info.coords.resize(3,2*tstruc.mol_info.names.size());


          {
            std::string struc_path=map_dir + "/parent_arrows.vasp";
            tstruc2.lat_column_mat=super_pstruc.lat_column_mat;
            tstruc2.mol_info.coords.leftCols(tstruc.mol_info.size())=super_pstruc.mol_info.coords;
            tstruc2.atom_info.coords.leftCols(tstruc.atom_info.size())=super_pstruc.atom_info.coords;

            Eigen::Matrix3d F=super_pstruc.lat_column_mat*mapped_cstruc.lat_column_mat.inverse();
            tstruc2.mol_info.coords.rightCols(tstruc.mol_info.size())=F*mapped_cstruc.mol_info.coords;
            tstruc2.atom_info.coords.rightCols(tstruc.atom_info.size())=F*mapped_cstruc.atom_info.coords;
            std::ofstream fout(struc_path);
            VaspIO::PrintPOSCAR POS(tstruc2);
            POS.print(fout);
          }              
          {
            std::string struc_path=map_dir + "/child_arrows.vasp";
            tstruc2.lat_column_mat=mapped_cstruc.lat_column_mat;
            tstruc2.mol_info.coords.leftCols(tstruc.mol_info.size())=mapped_cstruc.mol_info.coords;
            tstruc2.atom_info.coords.leftCols(tstruc.atom_info.size())=mapped_cstruc.atom_info.coords;

            Eigen::Matrix3d F=mapped_cstruc.lat_column_mat*super_pstruc.lat_column_mat.inverse();
            tstruc2.mol_info.coords.rightCols(tstruc.mol_info.size())=F*super_pstruc.mol_info.coords;
            tstruc2.atom_info.coords.rightCols(tstruc.atom_info.size())=F*super_pstruc.atom_info.coords;
            std::ofstream fout(struc_path);
            VaspIO::PrintPOSCAR POS(tstruc2);
            POS.print(fout);
          }
          {
            std::string struc_path=map_dir + "/final.vasp";
            std::ofstream fout(struc_path);
            VaspIO::PrintPOSCAR POS(super_cstruc);
            POS.print(fout);
          }
          {
            std::string struc_path=map_dir + "/as_mapped.vasp";
            std::ofstream fout(struc_path);
            VaspIO::PrintPOSCAR POS(mapped_cstruc);
            POS.print(fout);
          }
          if(interp_images){
            
            for(Index j=0; j<=interp_images; ++j){
            
              double cscale = double(j)/double(interp_images);
              double pscale=1.-cscale;
            
              std::string struc_path=map_dir + "/image." + std::to_string(j) + ".vasp";
            
              tstruc.lat_column_mat=pscale*super_pstruc.lat_column_mat+cscale*mapped_cstruc.lat_column_mat;
              tstruc.atom_info.coords=pscale*super_pstruc.atom_info.coords+cscale*mapped_cstruc.atom_info.coords;
              tstruc.mol_info.coords=pscale*super_pstruc.mol_info.coords+cscale*mapped_cstruc.mol_info.coords;
              std::ofstream fout(struc_path);
              VaspIO::PrintPOSCAR POS(tstruc);
              POS.print(fout);
            }
          }
        }
      }
    
    }
    ca_range[0]+=ca_range[1];
  }
  return 0;
}
