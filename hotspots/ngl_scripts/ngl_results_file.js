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


//create the NGL stage
document.addEventListener("DOMContentLoaded", function () {
    stage = new NGL.Stage("viewport", {backgroundColor: 'black'});
    load()
    get_slider_values()
    checkbox_toggles()
})