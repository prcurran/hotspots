/* Functions needed for running the NGL script */

function log(loaded_objects){
    console.log(loaded_objects)
    return loaded_objects
        }

function set_surface_isolevel(loaded_objects, grid_params){
    /*
    grid_params: {}
    grid_params[grid_type]: str, "donor", "acceptor", "apolar, "positive", "negative"
    grid_params[isolevel]: float, target contour level
    */
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

function set_dot_thresholds(loaded_objects, grid_params){
    /*
    grid_params: {}
    grid_params[grid_type]: str, "donor", "acceptor", "apolar", "positive", "negative"
    grid_params[threshold_max]: float, max contour level
    grid_params[threshold_min]: float, min contour level
    */
    for (var i=0; i< loaded_objects.length; i++){
        if (loaded_objects[i].name.includes(grid_params.grid_type)){
            for (var r=0; r<loaded_objects[i].reprList.length; r++){
                if(loaded_objects[i].reprList[r].name == "dot"){
                    loaded_objects[i].reprList[r].setParameters({thresholdMax: grid_params.threshold_max,
                                                                 thresholdMin: grid_params.threshold_min})
                }
            }
        }
    }

    for (var k=0;k<loaded_objects.length;k++){
        loaded_objects[k].updateRepresentations({position:true})
    }
    return loaded_objects
    }

function toggle_visibility(loaded_objects, params){
    /*params: {}
    params.visibility: bool, true for visible
    params.type: str, "surface", "dots", "cartoon", "lines", "liquorice", "ball and stick"
    params.selection: str, any other selections
    */

    if (typeof params.selection === "undefined"){
        params.selection=""}

    for (var i=0; i< loaded_objects.length; i++){
        if (loaded_objects[i].name.includes(params.selection)){
            for (var r=0; r<loaded_objects[i].reprList.length; r++){
                if(loaded_objects[i].reprList[r].name.includes(params.type)){
                    loaded_objects[i].reprList[r].setVisibility(params.visibility)
                    }
                }
            }
        }

    for (var k=0;k<loaded_objects.length;k++){
        loaded_objects[k].updateRepresentations({position:true})
        }
    return loaded_objects
    }


