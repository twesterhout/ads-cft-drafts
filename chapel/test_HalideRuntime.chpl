use HalideRuntime;

proc main() {

  var a = 1 .. 10;
  var b = 2 .. 11 by 3;
  writeln(a);
  writeln(b);
  writeln(rangeToDimension(a));
  writeln(rangeToDimension(b));
  writeln(dimensionToRange(new halide_dimension_t(0, 5, 2, 0)));

  var arr : [0 ..# 3] real;
  var buffer = new HalideBuffer(arr);

  var arr2 : [0 ..# 3, 0 ..# 2] real;
  var buffer2 = new HalideBuffer(arr2);
}
