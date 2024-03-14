import svgwrite
import re
import os

from .Graphical import *
from .Main_functions import *
from .Structure_Database import *

def check_pdb_path(pdb_file):
    """
    This method checks a path for a existing pdb-file with including a header (necessary for DSSP usage)

    Args:

    - pdb_file (str): Path to file that should be checked

    Returns:

    - the protein id from the file name if file is in correct form. Otherwise None
    
    """
    match = re.search(r'([a-zA-Z0-9]+)\.pdb$', pdb_file)
    #TODO check for header
    if match:
        base_name = os.path.basename(pdb_file)
        file_id = os.path.splitext(base_name)[0]
        with open(pdb_file, 'r') as input_file:
            pdb_content = input_file.read()
        if pdb_content.startswith("HEADER"):
            print("pdb-file with HEADER found!")
            return file_id
        else:
            print("pdb-file withput HEADER found")
            return None
    else:
        print(f"no pdb-file was found at: \"{pdb_file}\"")
        return None  

def visualize(dwg,protein,vis_type, AS_annotation, mark_endings, simple_helix, cysteins, simple_coil_avg, show_lddt_col,show_lddt_text):
    if len(protein.secondary_structures) == 0:
        return
    general_opacity=0.9
    only_path = False
    if(vis_type=='only-path'):
        #protein.get_protein_ordered_vis_objects(1, mark_endings)
        avg_path = create_simplified_path(protein.residues,averaging=simple_coil_avg)
        dwg.add(svgwrite.shapes.Polyline(points=avg_path, stroke='black', stroke_width=5, fill="none"))
        only_path=True
    
    elif vis_type=='simple-coil':
        protein.get_protein_ordered_vis_objects(simple_coil_avg, mark_endings)
        # vis non-coil ss noraml but connect by simplifying coil structure
        visualize_ordered_elements(dwg,protein.ordered_vis_elements, simple_helix, general_opacity, cysteins, show_lddt_col,show_lddt_text)

    elif vis_type=='normal':
        protein.get_protein_ordered_vis_objects(1, mark_endings)
        visualize_ordered_elements(dwg,protein.ordered_vis_elements, simple_helix, general_opacity=general_opacity,cysteins=cysteins, show_lddt_col=show_lddt_col,show_lddt_text=show_lddt_text)
    
    elif vis_type=='testing':
        protein.get_protein_ordered_vis_objects(1)
        create_testing_vis(dwg,ss_objects=protein.secondary_structures)
    
    elif vis_type=='fruchtermann':
        protein.get_protein_ordered_vis_objects(1)
        # vis of non coil ss segments that are pushe apart using Fruchtermann Reingold layout connected with simple lines
        do_fruchtermann_reingold_layout(protein.secondary_structures, k=0.5, iter=50)
        protein.scale_shift_coords(scale=1,shift=20, make_positive=True) #make everything positive
        visualize_ordered_elements(dwg,protein.ordered_vis_elements, simple_helix,general_opacity=general_opacity, cysteins=cysteins, show_lddt_col=show_lddt_col,show_lddt_text=show_lddt_text)
    
    else:
        raise ValueError("Error: Please provide a valid vis_type!")
        # additional:

    if AS_annotation:
        add_residue_annotations(dwg, protein.secondary_structures,only_path)

def get_best_rotation(pdb_file):
    """
    Calculates the best rotation for the input protein. The best rotation is based on maximising the area of content of the 2D representation, 
    minimising the number of overlapping non-coil secondary structures and minimising the depth (z range) of the picture.

    Args:

    - pdb_file (str): Path to pdb file the best rotation should be calculated for.

    Returns:

    - The rotation matrix that can be used to rotate the input protein in the best found orientation
    """
    file_id = check_pdb_path(pdb_file)
    protein = Protein()
    pdb_element = protein.parse_pdb(pdb_file)
    protein.get_secondary_structure(pdb_element,pdb_file)
    
    #rotate protein
    protein.find_best_view_angle()
    print()
    print("Best found rotation: ")
    print(protein.best_rotation)
    return protein.best_rotation

def db_get_SF_info(sf_number:int):
    """
    User can input a SCOP SF and get information on the SF in the database. 

    Args:

    - SF_numner(str): SCOP Superfamily identifier for that the information should be returned

    Returns:

    - The representative score, the protein-representative, the fixed rotation, and the fixed rotation type of the family
    """

    db = Structure_Database(None,None)
    return db.get_SF_info(int(sf_number))

def db_set_SF_pymol_rot(SF_number: str,pymol_output):
    """
    Can be used for setting a wanted fixed rotation for a specific Superfamily in the database. 
    Proteins matched to this SF will now be rotated using the new rotation matrix when using family_vis = True.
    The user can first rotate the family representative in pymol and then save the pymol rotation in the database. Every protein will than be rotated that way.
    
    Args:

    - SF_number (str): SCOP SuperFamily number that`s fixed rotation matrix should be changend.
    - pymol_output (str): Pymols get_view() output matrix. Contains rotation matrix
    
    Returns:

    - returns nothing but changes the db entry of the given SF_number. 
    """
    db= Structure_Database(None,None)
    rot_matrix = transform_pymol_out_to_UT(pymol_output)[0]
    db.set_manual_SF_rotation(SF_number, rot_matrix)

def create_2DSVG_from_pdb(pdb_file:str, result_dir:str, tmp_dir:str, family_vis:bool=True, fam_aligned_splitting:bool = True, drop_family_prob:float = 0.5,foldseek_executable:str=None,fixed_sf:int=None,fam_info:bool=True ,domain_splitted:bool=False, domain_annotation_file:str=None, domain_hull:bool=True, visualisation_type:str ="normal", 
                    cysteins:bool=True, as_annotation:bool=False, mark_endings:bool=True, simple_helix:bool=True, show_lddt_col:bool=False,show_lddt_text:bool=False, find_best_rot_step:int = 30, simple_coil_avg:int=10, chain_info:bool=True):
    """
    Main method of the package for creating 2D visualisations of a protein in form of a SVG file. The user can decide between different visualisation options.
    
    Args:\n\n

    - pdb_file (str): Path to pdb file the visualisation file is generated on. (Required)\n
    - result_dir (str): Path to dir, where the output file is saved (file name is automatically set based on input file name). (Required)\n
    - tmp_dir (str): Path to dir, where temporary files needed for analyis (e.g. foldseek) and visualisations are saved. (Required)\n 
    
    - family_vis (bool): If True, enables family-wise visualization, uses SCOP SF database with calculated representatives. Default is True.\n
    - fam_aligned_splitting (bool): If True, the protein is split into SF-aligned (is rotated based on this segment) and not-aligned parts. THey are connected with a dashed line. Default is True.\n
    - drop_family_prob (float): Allows the program to drop a chosen SF if the FoldSeek probability is smaller than given cut-off. In this case the protein rotation is determined using the implemented "find_best_view_angle" method. Default is 0.5. \n
    - show_lddt_col (bool): LDDT scores from FoldSeek alignment to best matching SF is shown per residue as colorscale (magenta). Default is False. \n
    - show_lddt_text (bool): LDDT scores from FoldSeek alignment to best matching SF is shown per residue. Default is False. \n
    - foldseek_executable (str): Path to foldseek executable (will be used for family alignment).\n
    - fixed_sf (int): Fixed SCOP superamily accession number that will be used for the protein (/ for every domain in the protein) (to distant SFs cannot be used for aligning)\n
    - fam_info (bool): If True, adds assigned Superfamily number and corresponding foldseek probability to the drawing. Default is True.\n

    - domain_splitted (bool): If True, protein is split into domains using the provided domain annotation file. Can be used in combination with family_vis which is then applied on each domain seperatly. Default is False.\n
    - domain_annotation_file (str): Path to the domain annotation file. Required if domain_splitted is True.\n
    - domain_hull (bool): If True sourounds domains with smoothed convex hull colored based on the secondary structure composition (R,G,B) <-> (helix,sheet,coil). Default is True\n

    - visualisation_type (string): "only-path", "normal", or "simple-coil". Default is "normal".\n
    - cysteins (bool): If True, includes calculated cystein bonds in the visualisation. Default is True.\n
    - as_annotation (bool): If True, includes AS-annotation. Default is False.\n
    - mark_endings (bool): If True, marks the endings. Default is True.\n
    - simple_helix (bool): If True, helix are represented in a simple way (file size eficient). Default is True.\n
    - find_best_rot_step (int): Is the size of steps per 3D rotation angle (total=3) taken to find the rotatoin showing the most of the input protein. Increasing leads to faster runtime but worse visualisations. Default is 30.\n
    - simple_coil_avg (int): Coil structures will be summarised together. e.g. 10 means that every 10 coil-residues will be averaged and treated as one point. Bigger values lead to higher simplification. Is only used when "simple-coil" or "only-path" is used. Default is 10\n
    - chain_info (bool): If true and multi-chain protein structure is given: adds chain annotations (from pdb) in the visualizations. Default is True.\n

    Returns: The path to the created SVG-file. \n

    - Creates a SVG file containing the 2D visualisation of the input protein in the given result_dir.
    """
    ############## Validate arguments ##############
    if domain_splitted and not domain_annotation_file:
        raise ValueError("Domain annotation file is required for domain-split analysis.")
    if visualisation_type=="split-alignment" and not family_vis:
        raise ValueError("Alignment split option can only be used when doing the family visualisation (aligned part is used for splitting).")
    valid_vis_types = ["only-path","normal","simple-coil"]
    if visualisation_type not in valid_vis_types:
        raise ValueError(f'"{visualisation_type}" is not a valid visualisation type. Please use one of the following: {valid_vis_types}')
    file_id = check_pdb_path(pdb_file)
    if file_id ==None:
        raise ValueError(f'"{pdb_file}" is not a valid pdb input. Please check the path and the file.')
    #counteract impossible combinations
    if not family_vis:
        fam_aligned_splitting=False
        drop_family_prob=False
    
    structure_database_obj = Structure_Database(foldseek_executable,tmp_dir)

    chain_pdb_dict = split_pdb_chains(pdb_file,tmp_dir)
    print(f"{len(chain_pdb_dict)} chains were found in the input pdb and will be visualized seperatly")
    print(chain_pdb_dict)
    chain_prots = []
    all_chain_domain_prots =[]
    for chain_id,chain_pdb in chain_pdb_dict.items():
        ############## 1) Split into domains if used ##############
        if domain_splitted:
            #TODO adapt domain splitting per chain...(chainsaw reformatting changes...)
            domain_files = get_domain_pdbs(chain_pdb,chain_id,domain_annotation_file,tmp_dir)
        else:
            domain_files = []
            domain_files.append(chain_pdb)

        print(f"{len(domain_files)} domain(s) were found for chain {chain_id} and will be visualized seperatly:")
        if len(domain_files) ==0:
            continue
        ############## 2) Get rotation of protein (family-vis / best rotation) ##############
        dom_proteins = []
        matched_sfs = []
        probs = []
        for dom_file in domain_files:
            print(f"\n### \"{dom_file}\" ###\n")
            if family_vis:
                sf,prob,sf_aligned_pdb_file, aligned_region,lddtfull = structure_database_obj.initial_and_fixed_Sf_rot_region(dom_file,drop_family_prob,fixed_sf=fixed_sf)
                if sf_aligned_pdb_file== None:
                    #no mathcing sf (with higher than min prob found) --> do normal vis
                    print("\033[91m"+ "No maching SF (with higher prob than min prob) found. Normal roation algo will be used!"+"\033[0m" )
                    add_header_to_pdb(dom_file)
                    dom_prot = Protein()
                    dom_prot.chain = chain_id
                    pdb_element = dom_prot.parse_pdb(dom_file)
                    dom_prot.get_secondary_structure(pdb_element,dom_file)
                    dom_prot.scale_shift_coords(scale=30,x_shift=0,y_shift=0, make_positive=False)
                    dom_prot.find_best_view_angle(step_width = find_best_rot_step)
                    dom_proteins.append(dom_prot)
                    matched_sfs.append("no_SF")
                    probs.append(-1)
                    continue
                    
                add_header_to_pdb(sf_aligned_pdb_file)
                dom_prot = Protein()
                dom_prot.chain = chain_id
                pdb_element = dom_prot.parse_pdb(sf_aligned_pdb_file)
                dom_prot.get_secondary_structure(pdb_element,sf_aligned_pdb_file)
                dom_prot.scale_shift_coords(scale=30,x_shift=0,y_shift=0, make_positive=False)
                dom_prot.add_lddt_to_aligned_residues(aligned_region,lddtfull) if show_lddt_col or show_lddt_text else None
                dom_prot.sf= sf
                matched_sfs.append(sf)
                dom_prot.prob= prob
                probs.append(prob)
                if fam_aligned_splitting:
                    # split protein in 3 segments (aligment-based)
                    front_part,aligned_part,end_part = dom_prot.split_aligned_part(aligned_region)
                    # shift left and right part to the sides and make positive again
                    front_aligned_x_shift = calc_x_overlap_distance_between_prots(front_part,aligned_part) + 200
                    aligned_end_x_shift = calc_x_overlap_distance_between_prots(aligned_part,end_part) + 200
                    end_part.scale_shift_coords(scale=1,x_shift=aligned_end_x_shift,y_shift=0,make_positive=False)
                    front_part.scale_shift_coords(scale=1,x_shift=-front_aligned_x_shift,y_shift=0,make_positive=False)
            else:
                # manually find best rotation and continue with that
                add_header_to_pdb(dom_file)
                dom_prot = Protein()
                dom_prot.chain = chain_id
                pdb_element = dom_prot.parse_pdb(dom_file)
                dom_prot.get_secondary_structure(pdb_element,dom_file)
                dom_prot.print_ss_objects()
                dom_prot.scale_shift_coords(scale=30,x_shift=0,y_shift=0, make_positive=False)
                dom_prot.find_best_view_angle(step_width = find_best_rot_step)
            dom_proteins.append(dom_prot)
        
            
        
        # push domains apart if necceray 
        #repeat_layout_pairwise(domains=dom_proteins, k=1, iter_steps=20)

        # make all coords postive
        chain_residues = [res for dom in dom_proteins for res in dom.residues]
        chain_prot = Protein()
        chain_prot.residues = chain_residues
        chain_prot.scale_shift_coords(scale=1,x_shift=0,y_shift=0,make_positive=True)
        #shift domains to be in linear line:
        shift_domains_in_x_line(dom_proteins, 100)
        #shift for border space
        #chain_prot.scale_shift_coords(scale=1,x_shift=100,y_shift=100,make_positive=False)
        chain_prots.append(chain_prot)
        all_chain_domain_prots.append(dom_proteins)
    
    full_residues = [res for chain_prot in chain_prots for res in chain_prot.residues]
    full_prot = Protein()
    full_prot.residues = full_residues
    #shift for border space
    full_prot.scale_shift_coords(scale=1,x_shift=100,y_shift=100,make_positive=False)
    shift_chains_in_y_line(chain_prots=chain_prots,chain_gap=140)
    
    ############## 3) & 4) Create visualisatoin as wanted: only-path / normal / simple-coil (+ fam_aligned_splitting) ##############
    viewbox = calculate_viewbox(full_prot.residues,300) 
    result_file_path = result_dir+"/"+file_id+'_'+visualisation_type+'_familyVis_'+str(family_vis)+'_simpleHelix_'+str(simple_helix)+'_vis.svg'
    dwg = svgwrite.Drawing(result_file_path, viewBox=viewbox)
    
    print("\n### Starting visualizing of protein (domains): ###")   
    for chain_doms in all_chain_domain_prots:
        last_dom=None
        for dom_prot in chain_doms:
            if fam_aligned_splitting and dom_prot.fam_aligned_parts!=None :
                front_part, aligned_part, end_part = dom_prot.fam_aligned_parts
                visualize(dwg,front_part,visualisation_type, as_annotation, mark_endings, simple_helix, cysteins, simple_coil_avg, show_lddt_col,show_lddt_text)
                visualize(dwg,aligned_part,visualisation_type, as_annotation, mark_endings, simple_helix, cysteins,simple_coil_avg, show_lddt_col,show_lddt_text)
                visualize(dwg,end_part,visualisation_type, as_annotation, mark_endings, simple_helix, cysteins,simple_coil_avg, show_lddt_col,show_lddt_text)
                add_dashed_line_between_proteins(dwg,front_part,aligned_part,end_part)
            else:
                visualize(dwg,dom_prot,visualisation_type, as_annotation, mark_endings, simple_helix, cysteins,simple_coil_avg, show_lddt_col,show_lddt_text)
            
            if domain_splitted:
                dom_prot.add_hull(dwg,dom_prot.get_hull_color(), opacity=0.25) if domain_hull else None
                dwg.add(last_dom.connect_to_protein_dashline(dom_prot)) if last_dom != None else None
            if fam_info:
                dom_prot.add_sf_info(dwg)
            last_dom=dom_prot
        if chain_info and len(all_chain_domain_prots)>1:
            chain_doms[0].add_chain_info(dwg)

            
    print(f"\n### Visualization done! SVG file was created at \"{result_file_path}\" ###")
    dwg.save()
    return result_file_path, matched_sfs,probs

