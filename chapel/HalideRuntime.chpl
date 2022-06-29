module HalideRuntime {

require "/home/tom/src/ads-cft/third_party/Halide-11.0.1-x86-64-linux/include/HalideRuntime.h";

// private use CPtr;
// private use SysCTypes;
private use CTypes;

extern record halide_dimension_t {
  var min : int(32);
  var extent : int(32);
  var stride : int(32);
  var flags : uint(32);
}

extern type halide_type_code_t = c_int;
extern const halide_type_int : halide_type_code_t;
extern const halide_type_uint : halide_type_code_t;
extern const halide_type_float : halide_type_code_t;
extern const halide_type_handle : halide_type_code_t;
extern const halide_type_bfloat : halide_type_code_t;

extern "struct halide_type_t" record halide_type_t {
  var code : uint(8);
  var bits : uint(8);
  var lanes : uint(16);
}

extern "struct halide_device_interface_t" record halide_device_interface_t {
  var device_malloc : c_fn_ptr;
  var device_free : c_fn_ptr;
  var device_sync : c_fn_ptr;
  var device_release : c_fn_ptr;
  var copy_to_host : c_fn_ptr;
  var copy_to_device : c_fn_ptr;
  var device_and_host_malloc : c_fn_ptr;
  var device_and_host_free : c_fn_ptr;
  var buffer_copy : c_fn_ptr;
  var device_crop : c_fn_ptr;
  var device_slice : c_fn_ptr;
  var device_release_crop : c_fn_ptr;
  var wrap_native : c_fn_ptr;
  var detach_native : c_fn_ptr;
  var compute_capability : c_fn_ptr;
  var impl : c_void_ptr;
}

extern record halide_buffer_t {
  /** A device-handle for e.g. GPU memory used to back this buffer. */
  var device : uint(64);
  /** The interface used to interpret the above handle. */
  var device_interface : c_ptr(halide_device_interface_t);
  /** A pointer to the start of the data in main memory. In terms of
   * the Halide coordinate system, this is the address of the min
   * coordinates (defined below). */
  var host : c_ptr(uint(8));
  /** flags with various meanings. */
  var flags : uint(64);
  /** The type of each buffer element. */
  var m_type : halide_type_t;
  /** The dimensionality of the buffer. */
  var dimensions : int(32);
  /** The shape of the buffer. Halide does not own this array - you
   * must manage the memory for it yourself. */
  var dim : c_ptr(halide_dimension_t);
  /** Pads the buffer up to a multiple of 8 bytes */
  var padding : c_void_ptr;
}


proc rangeToDimension(r : range(?t, BoundedRangeType.bounded, ?strided))
    where isIntegral(t) {
  return new halide_dimension_t(r.low:int(32), r.size:int(32), r.stride:int(32), 0);
}

proc dimensionToRange(d : halide_dimension_t) {
  return d.min ..# d.extent by d.stride;
}

proc toHalideType(type t) where t == real(32) {
  return new halide_type_t(halide_type_float:uint(8), 32, 1);
}

proc toHalideType(type t) where t == real(64) {
  return new halide_type_t(halide_type_float:uint(8), 64, 1);
}


record HalideBuffer {
  type eltType;
  param rank : int;
  var raw : halide_buffer_t;
  var _dims : c_array(halide_dimension_t, rank);

  proc init(arr : [] ?t) where arr.isRectangular() {
    this.eltType = t;
    this.rank = arr.rank;
    this.raw = new halide_buffer_t(
      device=0,
      device_interface=nil,
      host=c_ptrTo(arr[arr.domain.low]):c_ptr(uint(8)),
      flags=0,
      m_type=toHalideType(eltType),
      dimensions=rank,
      dim=nil,
      padding=nil
    );
    this.complete();
    for i in 0 ..# rank {
      this._dims[rank - 1 - i] = rangeToDimension(arr.dim(i));
    }
    this.raw.dim = _dims:c_ptr(halide_dimension_t);
  }
}


} // end module HalideRuntime
