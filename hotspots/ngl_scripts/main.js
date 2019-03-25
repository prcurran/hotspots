 /*Notes
    main.js:

    General: Each loaded file is put into a promise (https://developers.google.com/web/fundamentals/primers/promises), which is 
    stored in [loaded_promises] and propagated through final_promise.

    proteinRepresentation() and surfaceRepresentation() create an representation object (o), which holds the representations for the 
    protein or surface loaded. Objects are stored and propagated in [loaded_objects]. 

    loadSingleFile() forwards loaded file objects to protein/surfaceRepresentation() depending on the file extension.

    run() = sets up all structures necessary ([loaded_promises], [loaded_objects], final_promise). Also handles sliders (maybe change this later)

    runFunctionOnAllObjects() wrapper function needed to run the functions in custom_functions.js 
    func_to_run =function from custom_functions.js ; func_to_run updates [loaded_objects], runFunctionOnAllObjects updates final_promise
    

*/

//Functions needed for loading files into the viewer


function unzipBlob(blob,callback, loaded_promises, fname_list) {
  zip.useWebWorkers = false
  // use a zip.BlobReader object to read zipped data stored into blob variable
  zip.createReader(new zip.BlobReader(blob), function(zipReader) {
    // get entries from the zip file
    zipReader.getEntries(function(entries) {
      // get data from the first file
      console.log(entries)
      entries.forEach(function(element){
          //var fname = entries[i].filename
          //console.log(fname)
          element.getData(new zip.BlobWriter(), function(data){
                self=this;
                file=new File([data], self.fname)
                // close the reader and calls callback function with uncompressed data as parameter
                zipReader.close();
                callback(file,loaded_promises,fname_list)
          }.bind({fname: element.filename}))
          //console.log(a)
      });
     //zipReader.close();
    });
  }, onerror);
};

function onerror(message) {
  console.error(message);
}

    function loadSingleFile(f,loaded_promises, fname_list){
            console.log(f)
            var file_extension=f.name.split('.').pop()
            console.log(file_extension)
            //getBase64(file)
            if (file_extension=='pdb'){
                proteinRepresentation(f,loaded_promises)
            }else if (file_extension=='ccp4'){
                gridRepresentation(f,loaded_promises)
            }else if (file_extension == 'zip'){

                unzipBlob(f, loadSingleFile, loaded_promises, fname_list)
            }
            else{
                alert("Unknown file extension")
            }
       };


     function proteinRepresentation(f,loaded_promises){
            var promise_loaded;
            promise_loaded=stage.loadFile( f ).then(function (o) {
                o.addRepresentation("cartoon", { color: "lightgrey" })
                o.addRepresentation( "ball+stick", {
                    sele: "(( not polymer or hetero ) and not ( water or ion ))",
                    scale: 0.5
                } );
                o.autoView()
                return o;
              });

              loaded_promises.push(promise_loaded);

            };  

    function gridRepresentation(f,loaded_promises){
            if (f.name.includes("frequency")){
                var properties ={
                thresholdType:"value",
                thresholdMin:1,
                dotType:'sphere',
                color: 'black'
                }
                var promise_loaded
                promise_loaded=stage.loadFile( f ).then(function (o) {
                    o.addRepresentation("dot", properties);
                    return o;
                });
            loaded_promises.push(promise_loaded)
            return
            }

            else if(f.name.includes("ranges")){
                if (f.name.includes("acceptor")){
                    var colourscheme = "Reds"
                }else if (f.name.includes("donor")){
                    var colourscheme = "Blues"
                }else{
                    var colourscheme = "Purples"}

                var properties = {thresholdType:'value',
                thresholdMin: 0.000001,
                colorScheme:'value',
                colorScale:colourscheme,
                dotType:'sphere',
                radius: 0.3,
                opacity:1.0 };

                var promise_loaded
                promise_loaded=stage.loadFile( f ).then(function (o) {
                o.addRepresentation("dot", properties);
                //o.reprList[0].setVisibility(false)

                //o.autoView();

                return o;
                })
            loaded_promises.push(promise_loaded)
            return
            }

            else if(f.name.includes("donor")){
                var properties={color: "blue",
                isolevelType: "value",
                isolevel: 14.00,
                opacity: 0.9
                };
             }
             else if(f.name.includes("acceptor")){
                var properties={color: "red",
                isolevelType: "value",
                isolevel: 14.00,
                opacity: 0.9
                };
             }
             else if(f.name.includes("apolar")){
                var properties={color: "#FFF176",
                isolevelType: "value",
                isolevel: 14.00,
                //wireframe: true,
                //linewidth: 3,
                opacity: 0.2
                };
             }
             else{
                console.log('File is neither acceptor, donor or apolar grid')
                return
             }
            var promise_loaded;
            promise_loaded=stage.loadFile( f ).then(function (o) {
                o.addRepresentation("surface", properties);
                //o.autoView();
                return o;
             });

            loaded_promises.push(promise_loaded);

    };



      function get_slider_values(){

        var slider_values = {}

        document.getElementById("apolar_slider").oninput=function(){
            runFunctionOnAllObjects(set_isolevel,{"grid_type":"apolar", "isolevel": Number(this.value)})
            document.getElementById("apolar_label").innerHTML=this.value
            //console.log(this.value)
            slider_values.apolar = Number(this.value)
        }

        document.getElementById("acceptor_slider").oninput=function(){
            runFunctionOnAllObjects(set_isolevel,{"grid_type":"acceptor", "isolevel": Number(this.value)})
            document.getElementById("acceptor_label").innerHTML=this.value
            slider_values.acceptor = Number(this.value)
        }

        document.getElementById("donor_slider").oninput=function(){
            runFunctionOnAllObjects(set_isolevel,{"grid_type":"donor", "isolevel": Number(this.value)})
            document.getElementById("donor_label").innerHTML=this.value
            slider_values.donor = Number(this.value)
        }
        //console.log(slider_values)
        return slider_values
    }


    function run(fname_list){
    	window.addEventListener("change", function(){
            var all_files = $('#load-file')[0].files
            var loaded_promises=[]
            //for (var i=0; i< all_files.length; i++){
                //loadSingleFile(all_files[i],loaded_promises)
			//};

            if(typeof final_promise == "undefined"){

                for (var i=0; i< all_files.length; i++){
                    loadSingleFile(all_files[i],loaded_promises, fname_list)
                    fname_list.push(all_files[i].name)
                };
                //console.log(all_files.length)
                console.log(loaded_promises)
                console.log(fname_list)

                final_promise=Promise.all(loaded_promises).then(function(loaded_objects) {
                    console.log(loaded_objects);
                    return loaded_objects;
                });
                console.log(final_promise)

            } else {
                //in case further files are loaded
                console.log(all_files)
                console.log(loaded_promises)
                console.log("loading more objects")
                //console.log(fname_list)
                for (var i=0; i< all_files.length; i++){
                    var fname = all_files[i].name
                    console.log(typeof fname)
                    console.log(fname_list)
                    console.log(typeof fname_list[0])

                    if(!fname_list.includes(fname)){
                        console.log("file not in loaded")
                        loadSingleFile(all_files[i],loaded_promises, fname_list)
                        fname_list.push(fname)
                    }
                };
                loaded_promises.push(final_promise);
                final_promise=Promise.all(loaded_promises).then(function(loaded_objects) {
                    for (var i=0; i<loaded_objects.length; i++){
                        if(loaded_objects[i].constructor === Array){
                            var old_array = loaded_objects[i]
                            loaded_objects.splice(i, 1)
                            break;
                        }
                    }

                    for(var j=0; j<old_array.length; j++){
                        loaded_objects.push(old_array[j])
                    }
                    //console.log(loaded_objects.length)
                    //console.log(all_files.length)

                    return loaded_objects;
                });
            }
        });
        slider_values = get_slider_values()
        //console.log(slider_values)
    };




//create the NGL stage
    document.addEventListener("DOMContentLoaded", function () {
        stage = new NGL.Stage("viewport", {backgroundColor: 'white'});
        var fname_list = []
        run(fname_list);
        });

//Attempt to dynamically resize the NGL viewer:
window.onresize = function(){
    var viewport = document.getElementById('viewport');
    viewport.style.width = "100%";
    viewport.style.height = "100%";
};

   

//Function that propagates the final promise when loaded_objects is modified. Needed for run custom_functions

   function runFunctionOnAllObjects(func_to_run, external_parameter){
    //func_to_run - any function from custom_functions.js
    //object_array = [loaded_objects]
    try{
        final_promise=final_promise.then(function(object_array){
            if (external_parameter!=null){
                object_array=func_to_run(object_array,external_parameter)
            }else{
                object_array=func_to_run(object_array)
            }
            return object_array
        });
    }
    catch (err){
        console.log('You need to put a comma after the function name')
        return object_array
    }
   }

//Utility fxns used commonly:
function arrayMax(arr) {
  var len = arr.length, max = -Infinity;
  while (len--) {
    if (arr[len] > max) {
      max = arr[len];
    }
  }
  return max;
};




// Some utility fxns that are not really used


    function getBase64(f){
        var reader = new FileReader();
        reader.readAsDataURL(f);
        reader.onload = function () {
            console.log(reader.result);
        };
    }

    function dataURLtoFile(dataurl, filename) {
    var arr = dataurl.split(','), mime = arr[0].match(/:(.*?);/)[1],
        bstr = atob(arr[1]), n = bstr.length, u8arr = new Uint8Array(n);
    while(n--){
        u8arr[n] = bstr.charCodeAt(n);
    }
    return new File([u8arr], filename, {type:mime});
    }


    function download_blackjack(filename, text) {
        var pom = document.createElement('a');
        pom.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
        pom.setAttribute('download', filename);

        if (document.createEvent) {
            var event = document.createEvent('MouseEvents');
            event.initEvent('click', true, true);
            pom.dispatchEvent(event);
        }
        else {
            pom.click();
        }
    }


