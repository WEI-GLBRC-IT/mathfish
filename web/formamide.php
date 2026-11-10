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
                <li><a href="mismatch.html">Mismatch Analysis<span>&rarr;Evaluation of mismatch stability and its affect on formamide dissociation for given probe, perfect match target, and mismatched target triplet.</span></a></li>
                <li id ="selected"><a href="formamide.html">Formamide Curve Generator<span>&rarr;Derivation of formamide dissociation profiles for multiple probe/target pairs.</span></a></li>
                <li><a href="competitor.html">Competitor Analysis<span>&rarr;The effect of a competitor oligo on the formamide dissociation profile of the probe with both target and non-target.</span></a></li>
            </ul><ul id="vert2">
            </ul><ul id="vert2">
                <li><a href="helpdocuments.html">Help<span>&rarr; Guidelines for mathFISH</span></a></li>
                <li><a href="link.html">Links<span>&rarr; Links to useful websites</span></a></li>
                <li><a href="citation.html">Citations<span>&rarr;List of Citations</span></a></li>
                <li><a href="contact.php">Contact Us<span>&rarr; Give us your valuable feedback</span></a></li>
            </ul>
        </div>
     

<div id="content">

<?php
$var2len = [
    'organism' => 3500,
    'probe'    => 50,
];

foreach ($var2len as $var => $maxlen) {
    if (!empty($_POST[$var]) && is_array($_POST[$var])) {
        foreach ($_POST[$var] as $idx => $value) {
            if (!is_string($value) || strlen($value) > $maxlen) {
                echo "Long sequence error!! $var sequence is not allowed to be longer than $maxlen";
                exit;
            }
        }
    }
}
?>


            <H2> RESULTS : Formamide Curve Generator </H2>
<?php  
include 'PHPfunctions.php';
$dummy=shell_exec('python3 cleanup.py');

$varNames = ['T', 'NA','Po','targetmolecule', 'domain'];

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

// Parse selected rows: only indices 0-4, each at most once.
// WHY: formamide.html has 5 probe/target rows, and each index is one full ClustalW + UNAFold
// pass in formamide.py. Without this, one anonymous POST of "0 0 0 ..." queued unlimited
// passes and held a CPU and an Apache worker until they finished (SYS-4247).
$chkdata = $_POST['chkdata'] ?? '';
$checked_indices = [];
foreach (preg_split('/\s+/', trim(is_string($chkdata) ? $chkdata : '')) as $tok) {
    if (preg_match('/^[0-4]$/', $tok) && !in_array((int) $tok, $checked_indices, true)) {
        $checked_indices[] = (int) $tok;
    }
}
$organisms = is_array($_POST['organism'] ?? null) ? array_values($_POST['organism']) : [];
$probes = is_array($_POST['probe'] ?? null) ? array_values($_POST['probe']) : [];

// Map checked rows to the posted sequences.
// The browser disables unchecked rows, and disabled fields are not posted, so probe[] holds
// ONLY the checked rows, in row order: the k-th checked row is $probes[k]. Looking them up
// by row number instead dropped pairs whenever row 1 was unchecked or rows had gaps
// (rows 2+4 gave 1 pair, row 3 alone gave an error page). A client that posts every row
// still gets the row-number lookup.
if (count($probes) === count($checked_indices)) {
    $rows = array_keys($checked_indices);
} else {
    $rows = $checked_indices;
}

// Write probes
fwrite($file, "VAR_probes\n");
$selectedProbes = [];
foreach ($rows as $idx) {
    if (isset($probes[$idx])) {
        $selectedProbes[] = fstripCommandTricks($probes[$idx]);
    }
}
fwrite($file, implode(";", $selectedProbes) . "\n");

// Write organisms
fwrite($file, "VAR_organisms\n");
$selectedOrganisms = [];
foreach ($rows as $idx) {
    if (isset($organisms[$idx])) {
        $selectedOrganisms[] = fstripCommandTricks($organisms[$idx]);
    }
}
fwrite($file, implode(";", $selectedOrganisms) . "\n");

// Close the file
fclose($file);

// Run python
$cmd='python3 formamide.py "'.$randomString.'"';
$output=shell_exec($cmd);
//echo $cmd."<br>";
// echo $output;

// Collect the output
$filePath = "userFiles/outputVars_".$randomString. ".txt";

if (!file_exists($filePath) || filesize($filePath) == 0) {
    echo "<span style='font-size: 14px; color: red;'>Error encountered. Please email the developer with code: <b>$randomString</b></span>";
    exit;
}
$content = file_get_contents($filePath);
$vars = explode(";", $content);

// A short result means the script failed part-way; PHP 8 warns on missing keys.
if (count($vars) < 2) {
    echo "<span style='font-size: 14px; color: red;'>Error encountered. Please email the developer with code: <b>$randomString</b></span>";
    exit;
}

// echo $vars[0]."<br>";
// echo $vars[1]."<br>";

?>

            <br>
            <div id ="alignment">
				<?php
				// Define color list
				$colors = ['black', 'red', 'green', 'orange', 'magenta'];
				
				for ($i = 0; $i < min(5, count($vars) - 2); $i++) {
					$color = $colors[$i];
					echo '<font face="Courier New" size="4" style="color: ' . $color . ';">';
					echo '<b>PROBE ALIGNMENT WITH TARGET ORGANISM - ' . ($i + 1) . ':</b><br>';
					echo $vars[$i].'<br><br>';
					echo '</font>';
				}
				?>
			</div>
			
			

			<div id="warning" style="font-size: small; font-family: 'Courier New', Courier, monospace;">
                <div class="help">
                    <!-- <br> See <a href="helpdocuments.html#mismatch%20analysis" target="_blank">help</a> for the interpretation of results. -->
					<br><?php echo $vars[count($vars) - 2]; ?>
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

