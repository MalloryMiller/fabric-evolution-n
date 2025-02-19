let base_file = "../"

let target_frame = ".png"

let TEMP_RANGE = [-30.0, -2];
let TEMP_STEPS = 2.0;


let L = "L20"


let TAR_RANGES = {
    "cc": [-1, -0.000],
    "uc": [-1, -0.000],
    "ue": [0, 5],
    "ss": [0, 70],
}

let TAR_STEPS = {
    "cc": 0.2,
    "uc": 0.2,
    "ue": 1,
    "ss": 14,
}



function setup() {
    set_ranges()
}


function set_ranges() {

    let gambar = document.getElementById("gam");
    let tarbar = document.getElementById("tar");
    let gamtyp = document.getElementById("gam-t");
    let tartyp = document.getElementById("tar-t");

    let default_exp = document.getElementById("exp").value;
    gambar.max = TEMP_RANGE[1]
    gambar.min = TEMP_RANGE[0]
    gambar.value = TEMP_RANGE[0]
    gambar.step = TEMP_STEPS

    tarbar.max = TAR_RANGES[default_exp][1]
    tarbar.min = TAR_RANGES[default_exp][0]
    tarbar.value = TAR_RANGES[default_exp][0]
    tarbar.step = TAR_STEPS[default_exp]

    gamtyp.max = TEMP_RANGE[1]
    gamtyp.min = TEMP_RANGE[0]
    gamtyp.value = TEMP_RANGE[0]
    gamtyp.step = TEMP_STEPS

    tartyp.max = TAR_RANGES[default_exp][1]
    tartyp.min = TAR_RANGES[default_exp][0]
    tartyp.value = TAR_RANGES[default_exp][0]
    tartyp.step = TAR_STEPS[default_exp]

    console.log(default_exp)
    if (default_exp == "ss") {
        document.getElementById("target-units").innerHTML = "Angle (°)"
    } else {
        document.getElementById("target-units").innerHTML = "Change"

    }

    update()

}


function update_from_text() {
    document.getElementById("gam").value = document.getElementById("gam-t").value;
    document.getElementById("tar").value = document.getElementById("tar-t").value;

    update()
}



function update() {    
    document.getElementById("gam-t").value = document.getElementById("gam").value;
    document.getElementById("tar-t").value = document.getElementById("tar").value;
    settings = ""

    if (document.getElementById("corrected").checked) {
        settings += "solutions-vertical/";
        document.getElementById("corrected-label").setAttribute("class", "")
    } else {
        settings += "solutions/";
        document.getElementById("corrected-label").setAttribute("class", "disabled")
    }
    
    if (document.getElementById("ver").checked) {
        settings += "frames-isolated/";
        document.getElementById("ver-label").setAttribute("class", "")
    } else {
        settings += "frames/";
        document.getElementById("ver-label").setAttribute("class", "disabled")
    }



    settings += document.getElementById("exp").value
    if (document.getElementById("gam-on").checked) {
        settings += "/temp" + document.getElementById("gam").value;
        document.getElementById("gam").disabled = false;
        document.getElementById("gam-t").disabled = false;
    } else {
        document.getElementById("gam").disabled = true;
        document.getElementById("gam-t").disabled = true;
    }

    


    settings += "/lam/" + L + "/"
    settings += Number(document.getElementById("tar").value).toFixed(6);

    console.log(base_file + settings + target_frame);
    document.getElementById("chart").src = base_file + settings + target_frame;
}
