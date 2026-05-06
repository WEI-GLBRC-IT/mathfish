<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">

<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
        <title>Calculations</title>
        <link href="style.css" rel="stylesheet" type="text/css">
        <script language = "JavaScript">
            function warn(number)
            {
                if(number ==1)
                {
                    if(confirm("Mismatches with Target Organism are greater than 3. Do u still want to continue ???"))
                    {

                    }
                    else
                    {
                        //window.location="probeaff.html";
                        history.back();
                    }
                }
                else
                {
                    if(confirm("Mismatches with Non Target Organism are greater than 3. Do u still want to continue ???"))
                    {

                    }
                    else
                    {
                        //window.location="probeaff.html";
                        history.back();
                    }

                }

            }
        </script>

    </head>
    <body>
        <BR>
        <div id = "headimg">
            <img src="images/uwlogo_web_sm_ctr.png" height="100px" ALT="UW-Madison Logo">
        </div>
        <div class="ubercolortabs">
            <ul>
                <li><a href="index.html"><span>HOME</span></a></li>
                <li><a href="helpdocuments.html"><span>HELP</span></a></li>
                <li><a href="download.html"><span>DOWNLOAD</span></a></li>
                <li><a href="link.html"><span>LINKS</span></a></li>
                <li><a href="citation.html"><span>CITATIONS</span></a></li>
                <li><a href="contact.php"><span>CONTACT</span></a></li>
            </ul>
        </div>

        <BR>
        <br><br><div id="vert_menu">
            <!--
   Vertical menu with spans that contain hover description.
   Right-arrow is a separator for users not viewing style.
   Added blank li and span in beginning to fix alignment in IE.
            -->
            <ul id="vert"><li><span></span></li>
                <li><a href="probeaff.html"><strong>General Analysis</strong><span>&rarr;Complete  Analysis for any probe-target pair including reaction &#916;G<sup>o</sup>'s,  &#916;G<sup>o</sup>overall, hybridization efficiency, and formamide dissociation  profile.</span></a></li>
                <li><a href="dgovseries.html"><B>&#916;G<sup>o</sup><sub>overall</sub> Series</B><span>&rarr; Calculate &#916;G<sup>o</sup><sub>overall</sub>values for a series of probes of constant length targeting a given region on the organism of interest</span></a></li>
                <li><a href="deltaG1.html"><B>&#916;G<sup>o</sup><sub>1</sub> Series</B><span>&rarr; Calculate &#916;G<sup>o</sup><sub>1</sub>values for a series of probes of constant length targeting a given region on the organism of interest </span></a></li>
                <li><a href="deltaG2.html"><B>&#916;G<sup>o</sup><sub>2</sub> Series</B><span>&rarr; Calculate &#916;G<sup>o</sup><sub>2</sub>values for a series of probes of constant length targeting a given region on the organism of interest</span></a></li>
                <li><a href="deltaG3.html"><B>&#916;G<sup>o</sup><sub>3</sub> Series</B><span>&rarr; Calculate &#916;G<sup>o</sup><sub>3</sub>values for a series of probes of constant length targeting a given region on the organism of interest</span></a></li>
                <li id ="selected"><a href="mismatch.html">Mismatch Analysis<span>&rarr;Evaluation of mismatch stability and its affect on formamide dissociation for given probe, perfect match target, and mismatched target triplet.</span></a></li>
                <li><a href="formamide.html">Formamide Curve Generator<span>&rarr;Derivation of formamide dissociation profiles for multiple probe/target pairs.</span></a></li>
                <li><a href="competitor.html">Competitor Analysis<span>&rarr;The effect of a competitor oligo on the formamide dissociation profile of the probe with both target and non-target.</span></a></li>
            </ul><ul id="vert2">
            </ul><ul id="vert2">
                <li><a href="helpdocuments.html">Help<span>&rarr; Guidelines for mathFISH</span></a></li>
                <li><a href="link.html">Links<span>&rarr; Links to useful websites</span></a></li>
                <li><a href="citation.html">Citations<span>&rarr;List of Citations</span></a></li>
                <li><a href="contact.php">Contact Us<span>&rarr; Give us your valuable feedback</span></a></li>
            </ul>
        </div>

<?php // echo "102"; ?><?php // echo "130"; ?> 
<?php // echo "189"; ?>        

<div id="content">

<?php
$var2len = [
    'organism' => 3500,
	'mismatch' => 3500,
	'ipprobesequence' => 50,
];

foreach ($var2len as $var => $maxlen) {
    if (isset($_POST[$var]) && strlen($_POST[$var]) > $maxlen) {
        echo "Long sequence error!! $var is not allowed to be longer than $maxlen";
        exit;
    }
}
?>

            <H2> RESULTS : Mismatch Analysis </H2>
<?php 
// echo "193"; 
include 'PHPfunctions.php';
$dummy=shell_exec('python3 cleanup.py');
//echo "185"; 
$varNames = ['T', 'NA','Po','ipprobesequence','targetmolecule', 'domain','organism','mismatch'];

// Generate a random string for the file name
$randomString = generateRandomString();

// Construct the file path
$filePath = "userFiles/inputVars_" . $randomString . ".txt";

// Open the file for writing
$file = fopen($filePath, "w");

// Loop over variables
foreach ($varNames as $varName) {
	fwrite($file, "VAR_" . $varName . "\n");
    if (isset($_POST[$varName])) {		
        // Sanitize the input before writing to the file
        $sanitizedInput = fstripCommandTricks($_POST[$varName]);
		fwrite($file, $sanitizedInput . "\n");
    } else {
		fwrite($file, "NA\n");
	}
}

// Close the file
fclose($file);

// Run python
$cmd='python3 mismatch.py "'.$randomString.'"';
$output=shell_exec($cmd);
// echo $cmd."<br>";
// echo $output;

// Collect the output
$filePath = "userFiles/outputVars_".$randomString. ".txt";
$content = file_get_contents($filePath);

if (!file_exists($filePath) || filesize($filePath) == 0) {
    echo "<span style='font-size: 14px; color: red;'>Error encountered. Please email the developer with code: <b>$randomString</b></span>";
    exit;
}

$vars = explode(";", $content);

// echo $vars[0]."<br>";
// echo $vars[1]."<br>";

?>
            <br>
            <div id ="alignment">
                <font face= "Courier New" size ="4">
                    <b>PROBE ALIGNMENT WITH TARGET ORGANISM:</b><br>
<?php echo $vars[0]; ?>
                    <br><br><br>
                </font>


<font face= "Courier New" size ="4">
                    <b>PROBE ALIGNMENT WITH NON-TARGET ORGANISM:</b><br>
<?php echo $vars[8]; ?>
                    <br>
                </font>
					<br>
                    <br>
			</div>
			
			<div id="warning" style="clear: both;">
			    <H3>Thermodynamic Descriptors</H3>
                <table border="1" align="left" >
                    <tr>
                        <td></td>
                        <td><b>TARGET ORGANISM</b></td>
                        <td><b>NON-TARGET ORGANISM</b></td>
                        <td><b><div class="help"> &#916;Value<a href ="#">(?)<span>Values for Non-Target Organism - Values for Target Organism</span></a></div></b></td>
                    </tr>
                    <tr>
                        <td><B>&#916;G<sup>o</sup><sub>1</sub></B></td>
						<?php echo "<td>".$vars[2]."&nbsp;kcal/mol</td>"; ?>
						<?php echo "<td>".$vars[10]."&nbsp;kcal/mol</td>"; ?>
						<?php echo "<td>".$vars[16]."&nbsp;kcal/mol</td>"; ?>
                    </tr>
                    <tr>
                        <td><B>&#916;G<sup>o</sup><sub>2</sub></B></td>
						<?php echo "<td>".$vars[3]."&nbsp;kcal/mol</td>"; ?>
						<?php echo "<td>".$vars[11]."&nbsp;kcal/mol</td>"; ?>
						<td>NA</td>
                    </tr>
                    <tr>
                        <td><B>&#916;G<sup>o</sup><sub>3</sub></B></td>
						<?php echo "<td>".$vars[4]."&nbsp;kcal/mol</td>"; ?>
						<?php echo "<td>".$vars[12]."&nbsp;kcal/mol</td>"; ?>
						<?php echo "<td>".$vars[17]."&nbsp;kcal/mol</td>"; ?>
                    </tr>
                    <tr>
                        <td><B>&#916;G<sup>o</sup><sub>overall</sub></B></td>
						<?php echo "<td>".$vars[5]."&nbsp;kcal/mol</td>"; ?>
						<?php echo "<td>".$vars[13]."&nbsp;kcal/mol</td>"; ?>
						<?php echo "<td>".$vars[18]."&nbsp;kcal/mol</td>"; ?>
                    </tr>
                    <tr>
                        <td><B>[FA]<sub>m</sub>*</B></td>
						<?php echo "<td>".$vars[7]."&nbsp;%</td>"; ?>
						<?php echo "<td>".$vars[15]."&nbsp;%</td>"; ?>
						<?php echo "<td>".$vars[20]."&nbsp;%</td>"; ?>
                    </tr>
                    <tr>
                        <td><B>Hybridization Efficiency**</B></td>
						<?php echo "<td>".$vars[6]."</td>"; ?>
						<?php echo "<td>".$vars[14]."</td>"; ?>
						<?php echo "<td>".$vars[19]."</td>"; ?>
                    </tr>
                </table>
            </div>
			
			<div id="footnote" style="clear: both;padding-top: 8px;">
				* Melting formamide concentration<br>
				** At 0% formamide
			</div>

			<div id="warning" style="font-size: small; font-family: 'Courier New', Courier, monospace;">
                <div class="help">
                    <!-- <br> See <a href="helpdocuments.html#mismatch%20analysis" target="_blank">help</a> for the interpretation of results. -->
					<br><?php echo $vars[21]; ?>
                </div>
            </div>

		<div id="warning" style="clear: both;">
			<br>
            <H3>FORMAMIDE CURVE</H3>
			<?php
                $imagePath = "userFiles/faplot_" . $randomString . ".png";
                
                if (file_exists($imagePath)) {
                    echo '<img src="' . $imagePath . '" alt="Formamide Curve" />';
                } else {
                    echo '<p>Image not found!</p>';
                }
            ?>

            <br><br>


			 <!-- Button to download the already generated CSV file -->
			<form action="userFiles/facurve_<?php echo $randomString; ?>.csv" method="GET">
				<button type="submit">Download as CSV</button>
			</form>

            <br>
            <br>
        </div>

    </body>
</html>

