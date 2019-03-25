function savePDB(loaded_objects, params){
	//params = {fname:'', suffix:''}

	for (var i=0; i<loaded_objects.length; i++){
		if(loaded_objects[i].name.includes(params.fname)) {
			var newStructure = loaded_objects[i].structure
			console.log(newStructure)
		    var pdbWriter = new NGL.PdbWriter( newStructure, {});
		    var str_fname = loaded_objects[i].name.split('.')[0].concat(params.suffix);
		    pdbWriter.download(str_fname);
		    //console.log(pdbWriter.getData())
		}
	}
        return loaded_objects
}

function saveChain(loaded_objects, params){
    //params: fname = str
	//params: sele = chain name
	//params: suffix = str
	for (var i=0; i<loaded_objects.length; i++){
		if(loaded_objects[i].name.includes(params.fname)) {
			var newStructure = loaded_objects[i].structure
			var selection = new NGL.Selection(params.sele)

			//var str_fname = "balababulibo";
			var pdbWriter = new NGL.PdbWriter(newStructure, {});
			//pdbWriter.download(str_fname);
			pdbWriter._writeRecords()
			var re=new RegExp(/\s{1,}/)

			var newRecords=[];
			for (var m in pdbWriter._records){
				var splitted_line=pdbWriter._records[m].split(re)
				if (splitted_line[0]=="ATOM"){
					if (splitted_line[4].includes(params.sele)){
						newRecords.push(pdbWriter._records[m])
						continue
					}

				}//else if(splitted_line[0] =='HETATM'){
					//if(splitted_line[4].includes('A1613')){
					//newRecords.push(pdbWriter._records[i])
					//}
				//}
			}
			//console.log(pdbWriter._records.length)
			console.log(newRecords.length)
			//console.log(pdbWriter._records[10].split(re))
			var str_fname = loaded_objects[i].name.split('.')[0].concat(params.suffix);
			download_blackjack(str_fname.concat(".pdb"), newRecords.join('\n'))

		}
	}

	for (var k=0;k<loaded_objects.length;k++){
        loaded_objects[k].updateRepresentations({position:true})
    }
        return loaded_objects

}

function singleChain(chain_object){
	var newChain = chain_object.structure
	var str_fname = newChain.name.split('.')[0].concat("_A");
	var pdbWriter = new NGL.PdbWriter( newChain, {});
	//pdbWriter.download(str_fname);
	console.log(chain_object.atomCount)
	console.log(chain_object.structure)
}


function test_alignment(){
	var l, _i, _j, x, y, ali1, ali2
	ali1 = "-UTA21-BBA"
	ali2 = "MUT-21BBBA"
    i = 0
    j = 0
    n = 10
    var aliIdx1 = []
    var aliIdx2 = []
 
    for (l = 0; l < n; ++l) {
      x = ali1[ l ]
      y = ali2[ l ]
 
      _i = 0
      _j = 0
 
      if (x === '-') {
        aliIdx2[ j ] = false
      } else {
        aliIdx2[ j ] = true
        _i = 1
      }
 
      if (y === '-') {
        aliIdx1[ i ] = false
      } else {
        aliIdx1[ i ] = true
        _j = 1
      }
 
      i += _i
      j += _j
    }
 console.log(aliIdx1)
 console.log(aliIdx2)
}

function saveALignment(loaded_objects, known_bromo_list){
	//params = list of our structures (not dictionary)
    for (var i in loaded_objects){
        if(loaded_objects[i].name.includes(known_bromo_list[0]) && !loaded_objects[i].name.includes("ccp4")){
        	var struct1 = loaded_objects[i].structure;
        }
    }

    for (var i in loaded_objects){
        if(loaded_objects[i].name.includes(known_bromo_list[1]) && !loaded_objects[i].name.includes("ccp4")){
        	var struct2 = loaded_objects[i].structure;
        }
    }
     var _s1 = struct1.getView(new NGL.Selection(''));
     var _s2 = struct2.getView(new NGL.Selection(''));


    const seq1 = _s1.getSequence();
    const seq2 = _s2.getSequence();

    //gives "not a constructor error"
    var ali = new stage.Alignment(seq1.join(''), seq2.join(''));
    ali.calc();
    ali.trace();
    download_blackjack('align.fasta', ali);
}
