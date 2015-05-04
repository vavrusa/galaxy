var structs = {};

function getcolor(name) {
  var ccolor = null;
  var color_name = $('#color input:checked').attr('id');
  if (name != "main") {
      ccolor = color.uniform('red');
  } else {
    if (color_name == "color-uniform") {
      ccolor = color.uniform('grey');
    }
    if (color_name == "color-occupancy") {
      ccolor = color.byAtomProp('occupancy');
    }
    if (color_name == "color-tempfactor") {
      var color_map = ['red','white','blue']
      var threshold = parseFloat(1.0 - $('#color-gradient').val());
      var steps = [0.0, threshold, 1.0];
      ccolor = color.byAtomProp('tempFactor', color.gradient(color_map, steps));
    }
  }
  return ccolor;
}

function lines(name, struct) {
  viewer.rm(name)
  viewer.lines(name, struct, {color:  getcolor(name)});
}

function cartoon(name, struct) {
  viewer.rm(name)

  var ccolor = getcolor(name);
  if (ccolor == null) {
    ccolor = color.ssSuccession();
  }

  viewer.cartoon(name, struct, {color: ccolor});
}

function lineTrace(name, struct) {
  viewer.rm(name)
  viewer.lineTrace(name, struct, {color: getcolor(name)});
}

function sline(name, struct) {
  viewer.rm(name)
  viewer.sline(name, struct, {color: getcolor(name)});
}

function ballsAndSticks(name, struct) {
  viewer.rm(name)
  viewer.ballsAndSticks(name, struct, {color: getcolor(name)});
}

function trace(name, struct) {
  viewer.rm(name)
  viewer.trace(name, struct, {color: getcolor(name)});
}

function loadModel(file, name, func) {
  $.ajax({ url : file, success : function(data) {
    structs[name] = io.pdb(data);
    structure = structs[name];
    func(name, structs[name]);
    if (name == "main") {
        viewer.fitTo(structs[name]);
    }
    $('#sequence_box').show();
    if (file.indexOf('_limit') > -1) {
      $('#sequence_value').text('Sequence is not available for aggregated PDBs.');
    } else {
      $('#sequence_value').text(modelSequence(name));
    }
  }});
}

function modelSequence(name) {
    var seq = "";
    structs[name].eachResidue(function (res) {
      seq += res.name();
    });
    seq = seq.replace(/(.)(.)(.)/g, function (str, p1, p2, p3, off, s) {
      return p1.toUpperCase() + p2.toLowerCase() + p3.toLowerCase();
    });
    seq = seq.replace(/Ala/g, " A ");
    seq = seq.replace(/Asx/g, " B ");
    seq = seq.replace(/Cys/g, " C ");
    seq = seq.replace(/Asp/g, " D ");
    seq = seq.replace(/Glu/g, " E ");
    seq = seq.replace(/Phe/g, " F ");
    seq = seq.replace(/Gly/g, " G ");
    seq = seq.replace(/His/g, " H ");
    seq = seq.replace(/Ile/g, " I ");
    seq = seq.replace(/Lys/g, " K ");
    seq = seq.replace(/Leu/g, " L ");
    seq = seq.replace(/Met/g, " M ");
    seq = seq.replace(/Asn/g, " N ");
    seq = seq.replace(/Pro/g, " P ");
    seq = seq.replace(/Gln/g, " Q ");
    seq = seq.replace(/Arg/g, " R ");
    seq = seq.replace(/Ser/g, " S ");
    seq = seq.replace(/Thr/g, " T ");
    seq = seq.replace(/Val/g, " V ");
    seq = seq.replace(/Trp/g, " W ");
    seq = seq.replace(/Xaa/g, " X ");
    seq = seq.replace(/Tyr/g, " Y ");
    seq = seq.replace(/Glx/g, " Z ");
    seq = seq.replace(/\*\*\*/g, " * ");
    seq = seq.replace(/\s/g, "");
    return seq;
}