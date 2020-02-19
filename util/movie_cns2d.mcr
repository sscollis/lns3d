#!MC 800

# CNS2D movie maker

$!VarSet |MFBD| = './'
$!VarSet |data| = 'deriv.'
$!VarSet |start| = 000
$!VarSet |end|   = 4000 
$!VarSet |inc|   = 50
$!VarSet |NumFiles| = ( ( |end| - |start| ) / |inc| + 1 ) - 1

$!ExportSetup
  ExportFormat = Rastermetafile 
  ExportFName = "|MFBD|movie.rm"

$!DoubleBuffer On

$!Loop |NumFiles|
  $!DrawGraphics No
  $!VarSet |step| = ( |start| + ( |loop| - 1 ) * |inc| )
  $!ReadDataSet  '"|MFBD||data||step%.6d|.plt" "|MFBD|prim.|step%.6d|.plt"' ResetStyle = False
  $!DrawGraphics Yes
  $!RedrawAll
  $!DoubleBuffer Swap
  $!If |loop| == 1
    $!Export
      Append = No
  $!EndIf
  $!If |loop| > 1
    $!Export
      Append = Yes
  $!EndIf
$!EndLoop

$!DoubleBuffer off
