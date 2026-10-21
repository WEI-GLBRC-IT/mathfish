<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">

<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
        <title>Calculations</title>
        <link href="style.css" rel="stylesheet" type="text/css">
        <script language = "JavaScript">
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
                <li id ="selected"><a href="deltaG2.html"><B>&#916;G<sup>o</sup><sub>2</sub> Series</B><span>&rarr; Calculate &#916;G<sup>o</sup><sub>2</sub>values for a series of probes of constant length targeting a given region on the organism of interest</span></a></li>
                <li><a href="deltaG3.html"><B>&#916;G<sup>o</sup><sub>3</sub> Series</B><span>&rarr; Calculate &#916;G<sup>o</sup><sub>3</sub>values for a series of probes of constant length targeting a given region on the organism of interest</span></a></li>
                <li><a href="mismatch.html">Mismatch Analysis<span>&rarr;Evaluation of mismatch stability and its affect on formamide dissociation for given probe, perfect match target, and mismatched target triplet.</span></a></li>
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
if (!empty($_POST['organism']) && (!is_string($_POST['organism']) || strlen($_POST['organism']) > 3500)) {
    echo "Long sequence error!! organism is not allowed to be longer than 3500";
    exit;
}

if (!empty($_POST['pl']) && (int)$_POST['pl'] > 50) {
    echo "Long sequence error!! Probe length is not allowed to be greater than 50";
    exit;
}
?>

            <H2> RESULTS : &Delta;G<sup>o</sup><sub>2</sub> Series  </H2>
<?php 
// echo "193"; 
include 'PHPfunctions.php';
$dummy=shell_exec('python3 cleanup.py');
//echo "185"; 
$varNames = ['T', 'NA','organism','ts1','ts2','pl'];

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
$cmd='python3 dg2.py "'.$randomString.'"';
$output=shell_exec($cmd);
// echo $cmd."<br>";
// echo $output;

// Collect the output
$filePath = "userFiles/outputVars_".$randomString. ".txt";

if (!file_exists($filePath) || filesize($filePath) == 0) {
    echo "<span style='font-size: 14px; color: red;'>Error encountered. Please email the developer with code: <b>$randomString</b></span>";
    exit;
}
$content = file_get_contents($filePath);
$vars = explode(";", $content);

// Check if content contains "ERROR"
if (strpos($vars[0], "ERROR") !== false) {
    echo "<span style='font-size: 14px; color: red;'><b>$vars[0]</b></span>";
    exit;
}

// A short result means the script failed part-way; PHP 8 warns on missing keys.
if (count($vars) < 2) {
    echo "<span style='font-size: 14px; color: red;'>Error encountered. Please email the developer with code: <b>$randomString</b></span>";
    exit;
}

// echo $vars[0]."<br>";
// echo $vars[1]."<br>";

?>
            <div id ="warning">
                <?php echo $vars[0]; ?>
			</div>
			
			<div id="warning" style="font-size: small; font-family: 'Courier New', Courier, monospace;">
                <div class="help">
					<br><?php echo $vars[1]; ?>
                </div>
            </div>
			
		<div id="warning" style="clear: both;">
			<br>
            <H3>Walking Oligomer</H3>
			<?php
                $imagePath = "userFiles/dg2plot_" . $randomString . ".png";
                
                if (file_exists($imagePath)) {
                    echo '<img src="' . $imagePath . '" alt="dG2 series" />';
                } else {
                    echo '<p>Image not found!</p>';
                }
            ?>

            <br><br>


			 <!-- Button to download the already generated CSV file -->
			<form action="userFiles/dg2series_<?php echo $randomString; ?>.csv" method="GET">
				<button type="submit">Download as CSV</button>
			</form>

            <br>
            <br>
        </div>

    </body>
</html>

