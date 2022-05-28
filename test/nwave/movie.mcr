#!MC 800

# Primitive variable movie maker

$!VarSet |MFBD| = 'e:\codes\lns3d\test\nwave\'
$!VarSet |data| = 'deriv.'
$!VarSet |start| = 000
$!VarSet |end|   = 20 
$!VarSet |inc|   = 1
$!VarSet |NumFiles| = ( ( |end| - |start| ) / |inc| + 1 ) - 1

$!ExportSetup
  ExportFormat = Rastermetafile 
  ExportFName = "|MFBD|movie.rm"

$!DoubleBuffer On

$!Loop |NumFiles|
  $!DrawGraphics No
  $!VarSet |step| = ( |start| + ( |loop| - 1 ) * |inc| )
  $!If |step| < 10
    $!ReadDataSet '"|MFBD||data||step%.1d|.plt" "|MFBD|prim.|step%.1d|.plt"' ResetStyle = False
  $!Endif
  $!If |step| > 10 
    $!ReadDataSet '"|MFBD||data||step%.2d|.plt" "|MFBD|prim.|step%.2d|.plt"' ResetStyle = False
  $!EndIf
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
