<?php 
	function generateRandomString($length = 5) {
		$characters = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ';
		$charactersLength = strlen($characters);
		$randomString = '';
		for ($i = 0; $i < $length; $i++) {
			$randomString .= $characters[rand(0, $charactersLength - 1)];
		}
		return $randomString;
	}

	//ATTACK PREVENTION
	function fstripCommandTricks($string){
		$string=preg_replace('/"/','',$string);
		$string=preg_replace('/&/','',$string);
		$string=preg_replace('/\//','',$string);
		$string=preg_replace('/\|/','',$string);
		$string=preg_replace('/pkexec/','',$string);
		$string=preg_replace('/\>/','',$string);
		
		
		return $string;
	}
?>