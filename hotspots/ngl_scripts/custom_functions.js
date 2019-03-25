  

   KB = {
    "BAZ2BA": '((1893.CA or 1944.CA or 1888.CA or 1901.CA or 1887.CA) or (1893.N or 1944.N or 1888.N or 1901.N or 1887.N))',
    "4rvr" : '((1893.CA or 1944.CA or 1888.CA or 1901.CA or 1887.CA) or (1893.N or 1944.N or 1888.N or 1901.N or 1887.N))',
    "4xub": '((1893.CA or 1944.CA or 1888.CA or 1901.CA or 1887.CA) or (1893.N or 1944.N or 1888.N or 1901.N or 1887.N))',
    "BRD1A": '((53.CA or 54.CA or 59.CA or 67.CA or 110.CA) or (53.N or 54.N or 59.N or 67.N or 110.N))',
    "5enb": '((1344.CA or 1345.CA or 1350.CA or 1358.CA or 1401.CA) or (1344.N or 1345.N or 1350.N or 1358.N or 1401.N))',
    "5n49": '((585.CA or 586.CA or 591.CA or 599.CA or 642.CA) or (585.N or 586.N or 591.N or 599.N or 642.N))',
    "BRPF1": '((26.CA or 27.CA or 32.CA or 40.CA or 83.CA) or (26.N or 27.N or 32.N or 40.N or 83.N))',
    "4lc2" : '((702.CA or 703.CA or 708.CA or 716.CA or 759.CA) or (702.N or 703.N or 708.N or 716.N or 759.N))',

   };

   KB1 = {
   	 "BRD1A": '((53.CA or 54.CA or 59.CA or 67.CA or 110.CA) or (53.N or 54.N or 59.N or 67.N or 110.N))',
	 "BAZ2BA": '((1893.CA or 1944.CA or 1888.CA or 1901.CA or 1887.CA) or (1893.N or 1944.N or 1888.N or 1901.N or 1887.N))',
	 "4rvr" : '((1893.CA or 1944.CA or 1888.CA or 1901.CA or 1887.CA) or (1893.N or 1944.N or 1888.N or 1901.N or 1887.N))',
    	 "4xub": '((1893.CA or 1944.CA or 1888.CA or 1901.CA or 1887.CA) or (1893.N or 1944.N or 1888.N or 1901.N or 1887.N))',
   };

   KB2 = ['BRD1A', 'BAZ2BA']

   KB3 = {
    "BAZ2BA": '((1887.CA or 1888.CA or 1893.CA or 1901.CA or 1944.CA) or (1887.N or 1888.N or 1893.N or 1901.N or 1944.N))',
    "5enb": '((1339.CA or 1340.CA or 1345.CA or 1353.CA or 1396.CA) or (1339.N or 1340.N or 1345.N or 1353.N or 1396.N))'
   }
   
   KND1 = {
   	"DCP2B": "((118.CA or 143.CA or 147.CA or 148.CA or 193.CA) or (118.N or 143.N or 147.N or 148.N or 193.N))"
   };

    kinases = {
    'cdk2' : "((14.CA or 16.CA or 18.CA or 78.CA or 81.CA or 86.CA or 89.CA) or (14.N or 16.N or 18.N or 78.N or 81.N or 86.N or 89.N))",
    'gsk3b' : "((66.CA or 68.CA or 70.CA or 130.CA or 133.CA or 138.CA or 141.CA) or (66.N or 68.N or 70.N or 130.N or 133.N or 138.N or 141.N))"
    }


    function log(loaded_objects){
        console.log(loaded_objects)
        return loaded_objects
   }

   function paint_red(loaded_objects){
        loaded_objects[0].reprList[0].setColor("red")
        return loaded_objects
   }

   function set_isolevel(loaded_objects, grid_params){
   	//grid_params = {}
   	//grid_params[grid_type] = donor, acceptor, or apolar
   	//grid_params[isolevel] = isovalue you want to set it to
   	for (var i=0; i< loaded_objects.length; i++){
   		if (loaded_objects[i].name.includes(grid_params.grid_type)){
	   		for (var r=0; r<loaded_objects[i].reprList.length; r++){
	   			if(loaded_objects[i].reprList[r].name == "surface"){
	   				loaded_objects[i].reprList[r].setParameters({isolevel:grid_params.isolevel})
	   			}
	   		}
   		}
   	}

    for (var k=0;k<loaded_objects.length;k++){
        loaded_objects[k].updateRepresentations({position:true})
    }
   	return loaded_objects
   }

    function set_min_threshold(loaded_objects, grid_params){
   	//grid_params = {}
   	//grid_params[grid_type] = donor, acceptor, or apolar
   	//grid_params[threshold] = value you want to set it to
   	for (var i=0; i< loaded_objects.length; i++){
   		if (loaded_objects[i].name.includes(grid_params.grid_type)){
	   		for (var r=0; r<loaded_objects[i].reprList.length; r++){
	   			if(loaded_objects[i].reprList[r].name == "dot"){
	   				loaded_objects[i].reprList[r].setParameters({thresholdMin:grid_params.threshold})
	   			}
	   		}
   		}
   	}

    for (var k=0;k<loaded_objects.length;k++){
        loaded_objects[k].updateRepresentations({position:true})
    }
   	return loaded_objects
   }


   function superposeStructures(loaded_objects,known_bromo_list){
    //known_bromo_list: dictionary of the protein names[key] and the alignment selection[value]
    var structures_list={}
    for (var j in known_bromo_list){
        structures_list[j]=[]
        for (var i in loaded_objects){
            if(loaded_objects[i].name.includes(j) && !loaded_objects[i].name.includes("ccp4")){
                structures_list[j].push(loaded_objects[i].structure);
                console.log(structures_list)
            }
        }
    }

    var first_iter = true;
    for(var key in structures_list){
    	if(structures_list[key].length>0){
	        if(first_iter){
	        	var main_struct = structures_list[key][0];
	        	var main_selection = known_bromo_list[key];
	        	first_iter = false;
	        	continue;
	        }
            	NGL.superpose(structures_list[key][0], main_struct, false, known_bromo_list[key], main_selection);
        }
    }

    for(var m in structures_list){
    	if(structures_list[m].length>1){
    		var current_struct=structures_list[m][0];
    		var current_selection=known_bromo_list[m];
    		for(var r=1; r<structures_list[m].length; r++){
    			NGL.superpose(structures_list[m][r], current_struct, false, known_bromo_list[m], current_selection)
    		}
    	}
    }

    for (var k=0;k<loaded_objects.length;k++){
        loaded_objects[k].updateRepresentations({position:true})
    }
    return loaded_objects;
}


function colorAllStructures(loaded_objects, prot_and_col){
//prot_and_col = list of strings: first element= potein name, second = protein colour
//more accurate to say it colours non-surfaces
	for (var i in loaded_objects){
            if(loaded_objects[i].name.includes(prot_and_col[0]) && loaded_objects[i].reprList[0].name!='surface'){
            	for(var k in loaded_objects[i].reprList){
                	loaded_objects[i].reprList[k].setColor(prot_and_col[1]);
                	loaded_objects[i].reprList[k].setColor('element');
                }
            }
        }
    return loaded_objects;
}

function hideBallSticks(loaded_objects, prot_name){
//prot_name is optional (string); if given, will hide all ball+stick reprs of that protein
//otherwise, hides all balls and sticks reprs
	for (var i in loaded_objects){
            if(loaded_objects[i].name.includes(prot_name)){
            	for(var k in loaded_objects[i].reprList){
            		if(loaded_objects[i].reprList[k].name == 'ball+stick'){
                		loaded_objects[i].reprList[k].setVisibility(false);
                	}
                }
            }else if(typeof prot_name ==='undefined'){
            	for(var k in loaded_objects[i].reprList){
            		if(loaded_objects[i].reprList[k].name == 'ball+stick'){
                		loaded_objects[i].reprList[k].setVisibility(false);
                	}
                }
            }
        }
    return loaded_objects;
}

function hideCartoons(loaded_objects, prot_name){
//prot_name is optional (string); if given, will hide all ball+stick reprs of that protein
//otherwise, hides all balls and sticks reprs
	for (var i in loaded_objects){
            if(loaded_objects[i].name.includes(prot_name)){
            	for(var k in loaded_objects[i].reprList){
            		if(loaded_objects[i].reprList[k].name == 'cartoon'){
                		loaded_objects[i].reprList[k].setVisibility(false);
                	}
                }
            }else if(typeof prot_name ==='undefined'){
            	for(var k in loaded_objects[i].reprList){
            		if(loaded_objects[i].reprList[k].name == 'cartoon'){
                		loaded_objects[i].reprList[k].setVisibility(false);
                	}
                }
            }
        }
    return loaded_objects;
}

function hideLicorice(loaded_objects, prot_name){
//prot_name is optional (string); if given, will hide all ball+stick reprs of that protein
//otherwise, hides all balls and sticks reprs
	for (var i in loaded_objects){
            if(loaded_objects[i].name.includes(prot_name)){
            	for(var k in loaded_objects[i].reprList){
            		if(loaded_objects[i].reprList[k].name == 'licorice'){
                		loaded_objects[i].reprList[k].setVisibility(false);
                	}
                }
            }else if(typeof prot_name ==='undefined'){
            	for(var k in loaded_objects[i].reprList){
            		if(loaded_objects[i].reprList[k].name == 'licorice'){
                		loaded_objects[i].reprList[k].setVisibility(false);
                	}
                }
            }
        }
    return loaded_objects;
}


function hideLines(loaded_objects, prot_name){
//prot_name is optional (string); if given, will hide all ball+stick reprs of that protein
//otherwise, hides all balls and sticks reprs
	for (var i in loaded_objects){
            if(loaded_objects[i].name.includes(prot_name)){
            	for(var k in loaded_objects[i].reprList){
            		if(loaded_objects[i].reprList[k].name == 'line'){
                		loaded_objects[i].reprList[k].setVisibility(false);
                	}
                }
            }else if(typeof prot_name ==='undefined'){
            	for(var k in loaded_objects[i].reprList){
            		if(loaded_objects[i].reprList[k].name == 'line'){
                		loaded_objects[i].reprList[k].setVisibility(false);
                	}
                }
            }
        }
    return loaded_objects;
}

function showBinders(loaded_objects, params){
	// params =[protein_name, ligand_name(optional)]
	for (var i in loaded_objects){
            if(loaded_objects[i].name.includes(params[0]) && (typeof loaded_objects[i].structure!=='undefined')){
            	selection = "(( not polymer or hetero ) and not ( water or ion ))";
            	loaded_objects[i].addRepresentation('licorice', {sele: selection, colorScheme: "element", multipleBond: true});
            	loaded_objects[i].updateRepresentations({position:true});
               
            }
        }
    return loaded_objects;
}

function showAllResidues(loaded_objects, params){
	//params = ['protein_name', 'optional:selection']
	for (var i in loaded_objects){
            if(loaded_objects[i].name.includes(params[0]) && (typeof loaded_objects[i].structure!=='undefined')){
            	if (typeof params[1] !=='undefined'){
	            	selection = params[1];
	            	loaded_objects[i].addRepresentation('licorice', {sele: selection});
	            	continue;
            	}
            	selection = ' not(water or ion or hetero)';
            	loaded_objects[i].updateRepresentations({position:true});
               	loaded_objects[i].addRepresentation('licorice', {sele: selection, linewidth: 2 });
            }
        }
    return loaded_objects;
}


function isZero(i){
    return i == 0;
}

function showBindingResidues(loaded_objects, params) {
	for (var i in loaded_objects){
            if(loaded_objects[i].name.includes(params[0]) && (typeof loaded_objects[i].structure!=='undefined')){
            	var selection = new NGL.Selection("(( not polymer or hetero ) and not ( water or ion ))");
            	var radius = 5;
    			var atomSet = loaded_objects[i].structure.getAtomSetWithinSelection( selection, radius );
    			if (!atomSet._words.every(isZero)){
                    var atomSet2 =loaded_objects[i].structure.getAtomSetWithinGroup( atomSet );
                    var sele2 = atomSet2.toSeleString()+' and not(water or ion or hetero)';
                    loaded_objects[i].addRepresentation('licorice', {sele: sele2, colorScheme: 'element', multipleBond: true});
                    loaded_objects[i].updateRepresentations({position:true});
                }
               
            }
        }
    return loaded_objects;
}

function showBindingContacts(loaded_objects, params){
	//shows contacts withn 5A of the bound stuff
	for (var i in loaded_objects){
            if(loaded_objects[i].name.includes(params[0]) && (typeof loaded_objects[i].structure!=='undefined')){
            	var selection = new NGL.Selection("(( not polymer or hetero ) and not ( water or ion ))");
            	var radius = 5;
    			var atomSet = loaded_objects[i].structure.getAtomSetWithinSelection( selection, radius );
    			var atomSet2 =loaded_objects[i].structure.getAtomSetWithinGroup( atomSet );
    			var sele2 = atomSet2.toSeleString();
            	loaded_objects[i].addRepresentation('contact', {
            		masterModelIndex: 0,
                	weakHydrogenBond: true,
               		maxHbondDonPlaneAngle: 35,
               		linewidth: 1,
            		sele: sele2 + " or LIG"});   
            }
        }
    return loaded_objects;
}


 function superposeStuctNoSele(loaded_objects, known_bromo_list){
   // known_bromo_list = actually a list of names of the things
    var structures_list={}
    for (var j=0; j<known_bromo_list.length; j++){
        structures_list[known_bromo_list[j]]=[]
        for (var i in loaded_objects){
            if(loaded_objects[i].name.includes(known_bromo_list[j]) && !loaded_objects[i].name.includes("ccp4")){
                structures_list[known_bromo_list[j]].push(loaded_objects[i].structure);
            }
        }
    }
    console.log(structures_list)
    var first_iter = true;
    for(var key in structures_list){
        if(first_iter){
            var main_struct = structures_list[key][0];
            first_iter = false;
            continue;
        }
            NGL.superpose(structures_list[key][0], main_struct, true);   
    }

    for(var m in structures_list){
        if(structures_list[m].length>1){
            var current_struct=structures_list[m][0];
            for(var r=1; r<structures_list[m].length; r++){
                NGL.superpose(structures_list[m][r], current_struct, true);
                console.log("no selection specified")
                
            }
        }
    }

    for (var k=0;k<loaded_objects.length;k++){
        loaded_objects[k].updateRepresentations({position:true})
    }
    return loaded_objects;
}


function showWaters(loaded_objects, params){
    // params =[protein_name, color]
    for (var i in loaded_objects){
            if(loaded_objects[i].name.includes(params[0]) && loaded_objects[i].reprList[0].name!='surface'){
                selection = "water";
                if(params[1]!=="undefined"){
                    loaded_objects[i].addRepresentation('ball+stick', {sele: selection, color: params[1]});
                }
                else{
                    loaded_objects[i].addRepresentation('ball+stick', {sele: selection, colorScheme: "element"});
                }
                loaded_objects[i].updateRepresentations({position:true});
               
            }
        }
    return loaded_objects;
}

function showBindingSiteWaters(loaded_objects, params) {
    for (var i in loaded_objects){
            if(loaded_objects[i].name.includes(params[0]) && loaded_objects[i].reprList[0].name!='surface'){
                var selection = new NGL.Selection("(( not polymer or hetero ) and not ( water or ion ))");
                var radius = 10;
                var atomSet = loaded_objects[i].structure.getAtomSetWithinSelection( selection, radius );
                if (!atomSet._words.every(isZero)){
                    var sele2 = atomSet.toSeleString()+' and not (polymer or (( not polymer or hetero ) and not ( water or ion )))';

                    if(params[1]!=="undefined"){
                        loaded_objects[i].addRepresentation('ball+stick', {sele: sele2, color: params[1]});
                    }
                    else{
                        loaded_objects[i].addRepresentation('ball+stick', {sele: sele2, colorScheme: "element"});
                    }
                    loaded_objects[i].updateRepresentations({position:true});
                }
               
            }
        }
    return loaded_objects;
}

function showAllContacts(loaded_objects, params){
    //shows contacts withn 5A of the bound stuff
    //doesn't work very well
    for (var i in loaded_objects){
            if(loaded_objects[i].name.includes(params[0]) && (typeof loaded_objects[i].structure!=='undefined')){
                loaded_objects[i].addRepresentation('contact', {
                    masterModelIndex: 0,
                    weakHydrogenBond: true,
                    maxHbondDonPlaneAngle: 35,
                    linewidth: 1});   
            }
        }
    return loaded_objects;
}

function saveChainHetAtom(loaded_objects, params){
	//params: {fname: str,
	//         chain_sel : str,
	//         hetatom_sele : [],
	//         suffix: str}
	for (var i=0; i<loaded_objects.length; i++){
		if(loaded_objects[i].name.includes(params.fname)) {
		    var str_fname = loaded_objects[i].name.split('.')[0].concat(params.suffix);
			var newStructure = loaded_objects[i].structure
			var selection = new NGL.Selection(params.sele)

			//var str_fname = "balababulibo";
			var pdbWriter = new NGL.PdbWriter(newStructure, {});
			//pdbWriter.download(str_fname);
			pdbWriter._writeRecords()
			var re=new RegExp(/\s{1,}/)

			var newRecords=[];
			for (var i in pdbWriter._records){
				var splitted_line=pdbWriter._records[i].split(re)
				if (splitted_line[0]=="ATOM"){
					if (splitted_line[4].includes(params.chain_sel)){
						newRecords.push(pdbWriter._records[i])
						continue
					}

				}else if(splitted_line[0] =='HETATM'){
				    for (var n =0; n<params.hetatom_sele.length; n++){
				        if(splitted_line[4].includes(params.hetatom_sele[n])){
				        newRecords.push(pdbWriter._records[i])
				        }
				        else if (splitted_line[5] == params.hetatom_sele[n]){
				        newRecords.push(pdbWriter._records[i])
				        }
				    }
				}
			}
			console.log(pdbWriter._records.length)
			console.log(newRecords.length)
			//console.log(pdbWriter._records[10].split(re))
			download_blackjack(str_fname.concat(".pdb"), newRecords.join('\n'))

		}
	}
	return loaded_objects
}


function bitArrayUnion (arr1, arr2) {
    const words1 = arr1._words;
    const words2 = arr2._words;
    const count = Math.min(words1.length, words2.length);
    for (let k = 0; k < count; ++k) {
      words1[ k ] |= words2[ k ];
    }
    for (let k = words1.length; k < count; ++k) {
      words1[ k ] = 0;
    }
    return arr1;
  }


function bitArrayIntersect (arr1, arr2){
    const words1 = arr1._words;
    const words2 = arr2._words;
    //console.log(words1.length);
    //console.log(words2.length);
    const count = Math.min(words1.length, words2.length);
    for (let k = 0; k < count; ++k) {
      words1[ k ] &= words2[ k ];
    }
    for (let k = words1.length; k < count; ++k) {
      words1[ k ] = 0;
    }
    return arr1;
}


function makeint(str){
    return str = parseInt(str)}


function arrayIntersect(a, b)
//intersection function for regular arrays lifted from SO
{
  var ai=0, bi=0;
  var result = [];

  while( ai < a.length && bi < b.length )
  {
     if      (a[ai] < b[bi] ){ ai++; }
     else if (a[ai] > b[bi] ){ bi++; }
     else /* they're equal */
     {
       result.push(a[ai]);
       ai++;
       bi++;
     }
  }

  return result;
}


function fastSuperpose(loaded_objects, params){
// params.name = str(protein names)
// ONLY works if all PDBs have same res numbers
//figure a way to make it return residue numbers

    arr_list = [];
    for (var i in loaded_objects){
            if(typeof loaded_objects[i].structure!=='undefined'){
            	var selection = new NGL.Selection("(( not polymer or hetero ) and not ( water or ion ))");
            	var radius = 5;
    			var atomSet = loaded_objects[i].structure.getAtomSetWithinSelection( selection, radius );
    			if (!atomSet._words.every(isZero)){
                    var atomSet2 =loaded_objects[i].structure.getAtomSetWithinGroup( atomSet );
                    //Relies on backbone N being the first atom
                    var back_sel = new NGL.Selection("backbone and .N")
                    var backbone = loaded_objects[i].structure.getAtomSetWithinSelection(back_sel, 0);
                    console.log(atomSet2.toSeleString())
                    console.log(backbone.toSeleString())
                    atomSet2 = bitArrayIntersect(atomSet2, backbone);
                    console.log(atomSet2.toSeleString())
                    aa_index = [] ;
                    atoms = atomSet2.toSeleString().split("@")
                    atoms = atoms[1].split(",")
                    atoms = atoms.map(makeint)
                    //console.log(atoms)

                    for(let n=0; n<atoms.length; n++){

                        function isequal(val){
                            return val == atoms[n]
                            }

                        ind = loaded_objects[i].structure.residueStore.atomOffset.findIndex(isequal)

                        if (ind != -1){
                            aa_index.push(ind)
                            }
                        }

                    console.log(loaded_objects[i].name)
                    //console.log(aa_index)
                    //console.log(aa_index.length)
                    //console.log(atoms.length)

                    aa_list = [];

                    for (let u = 0; u< aa_index.length; u++){
                        aa_list.push(loaded_objects[i].structure.residueStore.resno[aa_index[u]])
                    }

                    console.log(aa_list)
                    //arr_list.push(atomSet2)
                    arr_list.push(aa_list)


                }
            }
        }

     for (var m=1; m < (arr_list.length); m++){
       arr_list[m] = arrayIntersect(arr_list[m-1], arr_list[m])
       console.log(arr_list[m])
       console.log(m)
     }

    //console.log(arr_list[m].toSeleString())

    //align_atoms = arr_list[m].toSeleString()
    console.log(arr_list[m-1])
    aa = arr_list[m-1]
    ca_list =[]
    n_list = []

    for (var i=0; i< aa.length; i++){
        ca_list.push(aa[i].toString() + ".CA")
        n_list.push(aa[i].toString()+".N")
    }
    sel_string = "((" + ca_list.toString() + ") or (" + n_list.toString() + "))"
    //console.log(sel_string)
    align_list = {[params.name] : sel_string}
    console.log(align_list)
    return loaded_objects = superposeStructures(loaded_objects, align_list);
    //return loaded_objects
}
