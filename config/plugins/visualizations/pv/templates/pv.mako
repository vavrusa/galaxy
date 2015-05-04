<%
    # Use root for resource loading.
    root = h.url_for( '/' )
%>
## ----------------------------------------------------------------------------

<!DOCTYPE HTML>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
${h.css( 'base', 'jquery-ui/smoothness/jquery-ui')}
${h.js( 'libs/jquery/jquery',
        'libs/jquery/jquery.migrate',
        'libs/jquery/jquery-ui')}
${h.stylesheet_link( root + 'plugins/visualizations/pv/static/pv.css' )}
${h.javascript_link( root + 'plugins/visualizations/pv/static/pv.min.js' )}
${h.javascript_link( root + 'plugins/visualizations/pv/static/pv-helpers.js' )}
</head>

## ----------------------------------------------------------------------------
<body>
    %if not embedded:
    <div id="chart-header">
        <h2 class="title">Protein Viewer on '${hda.name}' (entity "${hda.metadata.name}")</h2>
        <div id="tooltip"></div>
        <p class="subtitle">
            %if hda.metadata.chains:
            Molecule has ${hda.metadata.models} models, and following chains: ${hda.metadata.chains}.
            %endif
            <p>
                Information &mdash; <i>${hda.info}</i>
            </p>
        </p>
    </div>
    %endif

    <!-- UI elements -->
    <table style="width: 800px;">
      <tr>
        <td width="600">
          <div id="gl"></div>
        </td>
        <td>
          <div id="mol-style">
            <h3>Style</h3>
            <input type="radio" name="mol-style" id="style-cartoon" checked>Cartoon</input><br/>
            <input type="radio" name="mol-style" id="style-ballsandsticks">Balls &amp; Sticks</input><br/>
            <input type="radio" name="mol-style" id="style-lines">Lines</input><br/>
            <input type="radio" name="mol-style" id="style-line-trace">Line Trace</input><br/>
            <input type="radio" name="mol-style" id="style-sline">Smooth Line Trace</input><br/>
            <input type="radio" name="mol-style" id="style-trace">Trace</input><br/>
          </div>
          
          <div id="color">
            <h3>Color by</h3>
            <input type="radio" name="mol-color" id="color-uniform">Uniform</input><br/>
            <input type="radio" name="mol-color" id="color-succession">Succession</input><br/>
            <input type="radio" name="mol-color" id="color-occupancy">Occupancy</input><br/>
            <input type="radio" name="mol-color" id="color-tempfactor">Temperature factor</input><br/>
            <input type="range" id="color-gradient" step="0.01" min="0.01" max="0.99" />
          </div>

         <h3>Information</h3>
         <ul id="info">
           <li>Selected <span class="code" id="selection">(nothing)</span></li>
         </ul>
        </td>
      </tr>
      <tr>
        <td colspan="2" id="sequence_box" class="admonition admonition-notice">
        <p class="admonition-header" id="sequence_name">Molecule sequence</p>
        <textarea rows="1" class="admonition-content search-query" id="sequence_value">
        </textarea>
        </td>
      </tr>
    <tr>
      <td colspan="2" class="admonition admonition-tip">
      <p class="admonition-header">Note</p>
      <p class="admonition-content">
      You can rotate the model by dragging with the left mouse button and zoom by scrolling.
      You can display the residue number by left-clicking on the molecule.
      <p>
      </td>
    </tr>
    </table>

    <!-- Visualization code -->
    <script type="text/javascript">
    var hdaId   = '${trans.security.encode_id( hda.id )}',
        hdaExt  = '${hda.ext}',
        ajaxUrl = "${h.url_for( controller='/datasets', action='index')}/" + hdaId + "/display?to_ext=" + hdaExt;

    var viewer = pv.Viewer(document.getElementById('gl'), {
               width: 600,
               height: 600,
               antialias: true,
               quality: 'medium',
               outline: false
             });

    $('#style-cartoon').click(function() {
      cartoon("main", structs['main']);
    });
    
    $('#style-ballsandsticks').click(function() {
      ballsAndSticks("main", structs['main']);
    });
    
    $('#style-line-trace').click(function() {
      lineTracmaine("main", structs["main"]);
    });
    
    $('#style-lines').click(function() {
      lines("main", structs["main"]);
    });
    
    $('#style-trace').click(function() {
      trace("main", structs["main"]);
    });
    
    $('#style-sline').click(function() {
      sline("main", structs["main"]);
    });
    
    $('#style-tube').click(function() {
      tube("main", structs["main"]);
    });

    $('#color-succession').click(function() {
      $('#mol-style input:checked').click();
    });

    $('#color-tempfactor').click(function() {
      $('#mol-style input:checked').click();
    });
    
    $('#color-occupancy').click(function() {
      $('#mol-style input:checked').click();
    });

    $('#color-uniform').click(function() {
      $('#mol-style input:checked').click();
    });

    $('#color-gradient').on("change mousemove", function() {
      $('#mol-style input:checked').click(); 
      var estimate = $(this).val();
    });

    viewer.addListener("atomClicked", function(picked, originalEvent) {
      if (picked) {
        var atom = picked.object().atom;
        var res = atom.residue();
        var chain = res.chain();
        $('#selection').text(chain.name() + '/' + res.name() + res.num() + '/' + atom.name());
      }
    });

    viewer.addListener("atomDoubleClicked", function(picked, originalEvent) {
      if (picked === null) {
        viewer.fitTo(structs['main']);
        return;
      }
      var transformedPos = vec3.create();
      var newAtom = picked.object().atom;
      var pos = newAtom.pos();
      if (picked.transform()) {
          vec3.transformMat4(transformedPos, pos, picked.transform());
        viewer.setCenter(transformedPos, 0);
      } else {
        viewer.setCenter(pos, 0);
      }
      viewer.setZoom(50, 500);
    });

    $(document).ready(function() {
        $("#color input option:first").prop('checked', true);
        loadModel(ajaxUrl, 'main', cartoon);
    });
    </script>
</body>
