# LaTeX2HTML 2018 (Released Feb 1, 2018)
# Associate internals original text with physical files.


$key = q/Sec:para/;
$ref_files{$key} = "$dir".q|node6.html|; 
$noresave{$key} = "$nosave";

$key = q/SubSec:Examples/;
$ref_files{$key} = "$dir".q|node5.html|; 
$noresave{$key} = "$nosave";

$key = q/SubSec:para/;
$ref_files{$key} = "$dir".q|node7.html|; 
$noresave{$key} = "$nosave";

1;

