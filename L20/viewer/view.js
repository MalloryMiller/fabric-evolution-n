let base_file = "../solutions/frames/"

let target_frame = "50.png"

let GAM_RANGE = [1.0, 5.0];
let GAM_STEPS = 1.0;
let TAR_RANGE = [-0.95, 0.95];



let TAR_RANGES = {
    "cc": [-0.95, 0],
    "uc": [-0.95, 0.95],
    "ue": [1, 5],
    "ss": [10, 80],
}

let TAR_STEPS = {
    "cc": .1,
    "uc": .1,
    "ue": 1,
    "ss": 10,
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
    gambar.max = GAM_RANGE[1]
    gambar.min = GAM_RANGE[0]
    gambar.value = GAM_RANGE[0]
    gambar.step = GAM_STEPS

    tarbar.max = TAR_RANGES[default_exp][1]
    tarbar.min = TAR_RANGES[default_exp][0]
    tarbar.value = TAR_RANGES[default_exp][0]
    tarbar.step = TAR_STEPS[default_exp]

    gamtyp.max = GAM_RANGE[1]
    gamtyp.min = GAM_RANGE[0]
    gamtyp.value = GAM_RANGE[0]
    gamtyp.step = GAM_STEPS

    tartyp.max = TAR_RANGES[default_exp][1]
    tartyp.min = TAR_RANGES[default_exp][0]
    tartyp.value = TAR_RANGES[default_exp][0]
    tartyp.step = TAR_STEPS[default_exp]

    console.log(default_exp)
    if (default_exp == "ss") {
        document.getElementById("target-units").innerHTML = "Angle (°)"
    } else {
        document.getElementById("target-units").innerHTML = "Strain"

    }

    update()

}


function update_from_text() {
    console.log("updating from da text box")

    document.getElementById("gam").value = document.getElementById("gam-t").value;
    document.getElementById("tar").value = document.getElementById("tar-t").value;

    update()
}



function update() {
    console.log("UPDATING")
    
    document.getElementById("gam-t").value = document.getElementById("gam").value;
    document.getElementById("tar-t").value = document.getElementById("tar").value;


    settings = document.getElementById("exp").value
    if (document.getElementById("gam-on").checked) {
        settings += "/gam" + document.getElementById("gam").value;
        document.getElementById("gam").disabled = false;
        document.getElementById("gam-t").disabled = false;
    } else {
        document.getElementById("gam").disabled = true;
        document.getElementById("gam-t").disabled = true;
    }

    
    
    
    if (document.getElementById("lam").checked) {
        settings += "/lam0.15";
        document.getElementById("lamd-label").setAttribute("class", "")
    } else {
        document.getElementById("lamd-label").setAttribute("class", "disabled")
    }

    settings += "/target" + document.getElementById("tar").value + "/";
    //settings = "gam1.0_target0.26.nc/"

    console.log(base_file + settings + target_frame);
    document.getElementById("chart").src = base_file + settings + target_frame;
}
