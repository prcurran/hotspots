function get_slider_values(){

        var slider_values = {}


        document.getElementById("apolar_slider").oninput=function(){
            runFunctionOnAllObjects(set_surface_isolevel,{"grid_type":"apolar", "isolevel": Number(this.value)})
            runFunctionOnAllObjects(set_dot_threshold,{"grid_type":"apolar", "threshold_min": Number(this.value)})
            document.getElementById("apolar_label").innerHTML=this.value
            //console.log(this.value)
            slider_values.apolar = Number(this.value)
        }

        document.getElementById("acceptor_slider").oninput=function(){
            runFunctionOnAllObjects(set_surface_isolevel,{"grid_type":"acceptor", "isolevel": Number(this.value)})
            runFunctionOnAllObjects(set_dot_threshold,{"grid_type":"acceptor", "threshold_min": Number(this.value)})
            document.getElementById("acceptor_label").innerHTML=this.value
            slider_values.acceptor = Number(this.value)
        }

        document.getElementById("donor_slider").oninput=function(){
            runFunctionOnAllObjects(set_surface_isolevel,{"grid_type":"donor", "isolevel": Number(this.value)})
            runFunctionOnAllObjects(set_dot_threshold,{"grid_type":"donor", "threshold_min": Number(this.value)})
            document.getElementById("donor_label").innerHTML=this.value
            slider_values.donor = Number(this.value)
        }
        //console.log(slider_values)
        return slider_values
    }


    function checkbox_toggles(){
        document.getElementById("apolar_surface_checkbox").oninput=function(){
            runFunctionOnAllObjects(toggle_visibility, {visibility: this.checked, type:"surface", selection:"apolar"})
            }

        document.getElementById("acceptor_surface_checkbox").oninput=function(){
            runFunctionOnAllObjects(toggle_visibility, {visibility: this.checked, type:"surface", selection:"acceptor"})
            }

        document.getElementById("donor_surface_checkbox").oninput=function(){
            runFunctionOnAllObjects(toggle_visibility, {visibility: this.checked, type:"surface", selection:"donor"})
            }

        document.getElementById("apolar_dot_checkbox").oninput=function(){
            runFunctionOnAllObjects(toggle_visibility, {visibility: this.checked, type:"dot", selection:"apolar"})
            }

        document.getElementById("acceptor_dot_checkbox").oninput=function(){
            runFunctionOnAllObjects(toggle_visibility, {visibility: this.checked, type:"dot", selection:"acceptor"})
            }

        document.getElementById("donor_dot_checkbox").oninput=function(){
            runFunctionOnAllObjects(toggle_visibility, {visibility: this.checked, type:"dot", selection:"donor"})
            }
        }

    function ensurePromisesAreLoaded(parameter) {
        var loaded_promises=parameter;
        return new Promise(function (resolve, reject) {
            (function waitForPromises(){
                if (loaded_promises.length == expected_loaded_promises && loaded_promises.length > leftover_promises) return resolve();
                setTimeout(waitForPromises, 5);

                console.log("printing expected!!", expected_loaded_promises, loaded_promises.length)
                })();
            });
        }


    function load(){
    	document.getElementById("load-file").oninput=function(){
            var all_files = $('#load-file')[0].files
            var loaded_promises=[]
            //expected_loaded_promises= count_expected_loaded_promises(all_files)

            var load_new_files = new Promise(
				function(resolve,reject){
				    for (var i=0; i< all_files.length; i++){
				        loadSingleFile(all_files[i],loaded_promises)
						};
						resolve(loaded_promises.length);
				});


            if(typeof final_promise == "undefined"){
                expected_loaded_promises= count_expected_loaded_promises(all_files)
                leftover_promises = 0

				final_promise=load_new_files.then(function(fulfilled){
					final_promise=ensurePromisesAreLoaded(loaded_promises).then( function(fulfilled) {
						console.log(loaded_promises)
						final_promise=Promise.all(loaded_promises).then(function(loaded_objects) {
							console.log(loaded_objects);
							return loaded_objects;
						});
						return final_promise;
					});
					return final_promise;
				});

                console.log(final_promise)

            } else {
                //in case further files are loaded
                ``
                console.log("loading more objects")
                leftover_promises = 1
                expected_loaded_promises = leftover_promises + count_expected_loaded_promises(all_files)
                loaded_promises.push(final_promise);

                final_promise=load_new_files.then(function(fulfilled){
					final_promise=ensurePromisesAreLoaded(loaded_promises).then(function(fulfilled){

                        final_promise=Promise.all(loaded_promises).then(function(loaded_objects) {
                            console.log("in the looop!")
                            console.log(loaded_objects)
                            for (var i=0; i<loaded_objects.length; i++){
                                if(loaded_objects[i].constructor === Array){
                                    console.log(loaded_objects[i])
                                    var old_array = loaded_objects[i]
                                    loaded_objects.splice(i, 1)
                                    break;
                                    }
                                }
                            console.log(old_array)
                            for(var j=0; j<old_array.length; j++){
                                loaded_objects.push(old_array[j])
                                }

                            return loaded_objects;
                            })

                        return final_promise;
                        })
                    return final_promise;
                    })
                //*/

            }
        }
    }


//create the NGL stage
    document.addEventListener("DOMContentLoaded", function () {
        stage = new NGL.Stage("viewport", {backgroundColor: 'black'});
        var fname_list = []
        load();
        slider_values = get_slider_values()
        checkbox_toggles()
        });