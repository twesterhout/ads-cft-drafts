(import (chezscheme))

(load-shared-object "libc.so.6")
(define c-sleep
  (foreign-procedure __collect_safe "sleep" (unsigned) unsigned))

(load-shared-object "/home/tom/src/petsc-prefix/lib/libpetsc.so.3.16")

(define-ftype char* (* char))
(define-ftype char** (* char*))

(define PetscInitialize
  (let ([argc 0]
        [argv '()]
        [c-PetscInitialize (foreign-procedure "PetscInitialize"
                             ((* int) (* char**) string string) int)])
    c-PetscInitialize))

(define argv (command-line-arguments))
(display argv)
(newline)

(load-shared-object "/home/tom/src/ads-cft/kernels/build/libmy_runtime.so")
(load-shared-object "./libutils.so")

(define halide-malloc
  (foreign-procedure "custom_alloc" (uptr size_t) uptr))

(define halide-free
  (foreign-procedure "custom_free" (uptr uptr) void))

(define-ftype halide_type_t
  (struct
    [code unsigned-8]
    [bits unsigned-8]
    [lanes unsigned-16]))

(define print-halide-type
  (foreign-procedure "print_halide_type_t" ((* halide_type_t)) void))


(define-ftype halide_dimension_t
  (struct
    [min integer-32]
    [extent integer-32]
    [stride integer-32]
    [flags unsigned-32]))

(define-ftype halide_buffer_t
  (struct
    [device unsigned-64]
    [device_interface void*]
    ; host is actually uint8_t*
    [host uptr]
    [flags unsigned-64]
    [type halide_type_t]
    [dimensions integer-32]
    [dim (* halide_dimension_t)]
    [padding void*]))

(define-record-type tensor
  (fields storage dtype shape strides offset))

(define-ftype Storage
  (struct
    [refcount uptr]
    [data uptr]))
(define Storage-guardian (ftype-guardian Storage))
(define free-dropped-Storage
  (lambda ()
    (let ([s (Storage-guardian)])
      (when s
        (printf "freeing ~s\n" (ftype-ref Storage (data) s))
        (halide-free 0 (ftype-ref Storage (data) s))
        (foreign-free (ftype-pointer-address s))
        (free-dropped-Storage)))))
(define malloc-Storage
  (lambda (dtype n)
    (free-dropped-Storage)
    (let ([s (make-ftype-pointer Storage
                                 (foreign-alloc (ftype-sizeof Storage)))])
      (ftype-set! Storage (refcount) s 1)
      (ftype-set! Storage (data) s (halide-malloc 0 (* (foreign-sizeof dtype) n)))
      (printf "allocated ~s\n" (ftype-ref Storage (data) s))
      (Storage-guardian s)
      s)))
(define Storage-data
  (lambda (s)
    (ftype-ref Storage (data) s)))

(define product
  (lambda (xs) (fold-left * 1 xs)))

(define vector-drop
  (lambda (k xs)
    (let ([n (max 0 (- (vector-length xs) k))])
      (do ((ys (make-vector n))
           (i 0 (+ i 1)))
          ((>= i n) ys)
        (vector-set! ys i (vector-ref xs (+ i k)))))))

(define vector-reverse
  (lambda (xs)
    (do ((n (vector-length xs))
         (ys (make-vector (vector-length xs)))
         (i 0 (+ i 1)))
        ((= i n) ys)
      (vector-set! ys i (vector-ref xs (- (- n 1) i))))))

(define row-major-strides
  (lambda (shape)
    (list->vector
      (let loop ([acc (list 1)]
                 [k 1]
                 [ls (reverse shape)])
        (if (null? ls)
          (cdr acc)
          (let ([new_k (* k (car ls))])
            (loop (cons new_k acc) new_k (cdr ls))))))))

(define row-major?
  (lambda (shape strides)
    (let ([rank (vector-length shape)])
      (let loop ([sᵢ₋₁ (vector-ref strides 0)]
                 [i 1])
        (if (< i rank)
            (let ([n (vector-ref shape i)]
                  [sᵢ (vector-ref strides i)])
              (if (= (div sᵢ₋₁ sᵢ) n)
                  (loop sᵢ (+ i 1))
                  #f))
            #t)))))

(define tensor-row-major?
  (lambda (t)
    (row-major? (tensor-shape t) (tensor-strides t))))

(define tensor-new
  (lambda (dtype shape)
    (cond
      [(eq? dtype 'double-float)
         (let* ([size (product shape)]
                [storage (malloc-Storage dtype size)]
                [strides (row-major-strides shape)])
           (make-tensor storage dtype (list->vector shape) strides 0))])))

(define tensor-ndim
  (lambda (t) (vector-length (tensor-shape t))))
(define tensor-data
  (lambda (t) (+ (Storage-data (tensor-storage t)) (tensor-offset t))))

(define-syntax ftype-alloc
  (syntax-rules ()
   [(ftype-alloc dtype count)
    (make-ftype-pointer dtype (foreign-alloc (* (ftype-sizeof dtype) count)))]))

(define symbol->halide_type_t
  (lambda (symbol ptr)
    (cond
      [(eq? symbol 'double-float)
         (ftype-set! halide_type_t (code) ptr 2)
         (ftype-set! halide_type_t (bits) ptr 64)
         (ftype-set! halide_type_t (lanes) ptr 1)])))

(define tensor->halide_buffer_t
  (lambda (t)
    (let ([p (ftype-alloc halide_buffer_t 1)]
          [dims (ftype-alloc halide_dimension_t (tensor-ndim t))])
      (ftype-set! halide_buffer_t (device) p 0)
      (ftype-set! halide_buffer_t (device_interface) p 0)
      (ftype-set! halide_buffer_t (host) p (tensor-data t))
      (ftype-set! halide_buffer_t (flags) p 0)
      (symbol->halide_type_t (tensor-dtype t) (ftype-&ref halide_buffer_t (type) p))
      (ftype-set! halide_buffer_t (dimensions) p (tensor-ndim t))
      (ftype-set! halide_buffer_t (dim) p dims)
      (do ([n (tensor-ndim t)]
           [sizes (tensor-shape t)]
           [strides (tensor-strides t)]
           [i 0 (+ i 1)])
          ((= i n))
        (ftype-set! halide_dimension_t (min) dims i 0)
        (ftype-set! halide_dimension_t (extent) dims i (vector-ref sizes (- (- n 1) i)))
        (ftype-set! halide_dimension_t (stride) dims i (vector-ref strides (- (- n 1) i)))
        (ftype-set! halide_dimension_t (flags) dims i 0))
      (ftype-set! halide_buffer_t (padding) p 0)
      p)))

(define halide_buffer_t-free
  (lambda (b)
    (foreign-free (ftype-pointer-address (ftype-ref halide_buffer_t (dim) b)))
    (foreign-free (ftype-pointer-address b))))

(define tensor-free
  (lambda (x) (halide-free 0 (tensor-data x))))

(define with
  (lambda (acquire release body)
    (let* ([resource (acquire)]
           [result (body resource)])
      (release resource)
      result)))

(define alloca
  (lambda (nbytes action)
    (with (delay (foreign-alloc nbytes)) foreign-free action)))

(define-syntax with-ftype
  (syntax-rules ()
    [(with-ftype dtype action)
       (alloca
         (ftype-sizeof dtype)
         (lambda (p) (action (make-ftype-pointer dtype p))))]))

(define with-halide-double
  (lambda (action)
    (with-ftype halide_type_t
      (lambda (p)
        (ftype-set! halide_type_t (code) p 2)
        (ftype-set! halide_type_t (bits) p 64)
        (ftype-set! halide_type_t (lanes) p 1)
        (action p)))))

(trace-define-syntax pointer-get
  (syntax-rules ()
    [(_ ptr dtype i) (foreign-ref dtype ptr (* i (foreign-sizeof dtype)))]))

(define tensor-numel
  (lambda (t)
    (product (vector->list (tensor-shape t)))))

(define tensor-stride
  (lambda (t i)
    (vector-ref (tensor-strides t) i)))

(define tensor-size
  (lambda (t i)
    (vector-ref (tensor-shape t) i)))

(define tensor-get
  (case-lambda
    ((t i) (pointer-get (tensor-data t)
                        (tensor-dtype t)
                        (* i (tensor-stride t 0))))
    ((t i j) (pointer-get (tensor-data t)
                          (tensor-dtype t)
                          (+ (* i (tensor-stride t 0))
                             (* j (tensor-stride t 1)))))))

(define tensor-slice
  (lambda (t dim i)
    (let* ([dtype (tensor-dtype t)]
           [offset (+ (tensor-offset t)
                      (* i (tensor-stride t dim) (foreign-sizeof dtype)))]
           [ndim (- (tensor-ndim t) 1)]
           [shape (make-vector ndim)]
           [strides (make-vector ndim)])
      (do ([i 0 (+ i 1)])
          ((= i dim))
        (vector-set! shape i (tensor-size t i))
        (vector-set! strides i (tensor-stride t i)))
      (do ((i dim (+ i 1)))
          ((= i ndim) (make-tensor (tensor-storage t) dtype shape strides offset))
        (vector-set! shape i (tensor-size t (+ i 1)))
        (vector-set! strides i (tensor-stride t (+ i 1)))))))

(define tensor-reshape
  (lambda (t shape)
    (unless (= (product shape) (tensor-numel t))
      (assertion-violation 'tensor-reshape
                           "incompatible shapes"
                           shape
                           (tensor-shape t)))
    (cond
      [(tensor-row-major? t)
        (let* ([s (tensor-stride t (- (tensor-ndim t) 1))]
               [strides (vector-map (lambda (x) (* x s)) (row-major-strides shape))])
          (make-tensor (tensor-storage t)
                       (tensor-dtype t)
                       (list->vector shape)
                       strides
                       (tensor-offset t)))]
      [else (assertion-violation 'tensor-shape "only row-major tensors may be reshaped")])))

(define tensor->vector
  (lambda (t)
    (cond
      [(= (tensor-ndim t) 0) (list)]
      [(= (tensor-ndim t) 1)
        (let ([v (make-vector (tensor-size t 0))])
          (do ([i 0 (+ i 1)])
              ((= i (vector-length v)) v)
            (vector-set! v i (tensor-get t i))))]
      [else
        (let ([v (make-vector (tensor-size t 0))])
          (do ([i 0 (+ i 1)])
              ((= i (vector-length v)) v)
            (vector-set! v i (tensor->list (tensor-slice t 0 i)))))])))

(define differentiation_matrix_bounded_kernel
  (lambda (lower upper t)
    (let ([foreign-kernel (foreign-procedure "differentiation_matrix_bounded"
                            (double-float double-float int (* halide_buffer_t)) void)]
          [n (tensor-size t 0)])
      (with (delay (tensor->halide_buffer_t t)) halide_buffer_t-free
        (lambda (buffer)
          (foreign-kernel lower upper n buffer))))))

(define ∂-matrix-bounded
  (lambda (lower upper n)
    (let ([t (tensor-new 'double-float (list n n))])
      (differentiation_matrix_bounded_kernel lower upper t)
      t)))

(define grid-bounded
  (lambda (lower upper n)
    (let ([foreign-kernel (foreign-procedure "bounded_grid"
                            (double-float double-float int (* halide_buffer_t)) void)]
          [t (tensor-new 'double-float (list n))])
      (with (delay (tensor->halide_buffer_t t)) halide_buffer_t-free
        (lambda (buffer)
          (foreign-kernel lower upper n buffer)
          t)))))

(define-ftype
  [tblis-len_type ptrdiff_t]
  [tblis-stride_type ptrdiff_t]
  [tblis-label_type char]
  [tblis-type int]
  [tblis-scalar
    (struct
      (type tblis-type)
      (scalar (union (float single-float)
                     (double double-float)
                     (complex-float (array 2 single-float))
                     (complex-double (array 2 double-float)))))]
  [tblis-tensor
    (struct
      (type tblis-type)
      (conj int)
      (scalar tblis-scalar)
      (data void*)
      (ndim unsigned)
      (len (* tblis-len_type))
      (stride (* tblis-stride_type)))])

; typedef enum
; {
;     TYPE_SINGLE   = 0,
;     TYPE_FLOAT    = TYPE_SINGLE,
;     TYPE_DOUBLE   = 1,
;     TYPE_SCOMPLEX = 2,
;     TYPE_DCOMPLEX = 3
; } type_t;
(define symbol->tblis-type
  (lambda (dtype)
    (cond
      [(equal? dtype 'single-float) 0]
      [(equal? dtype 'double-float) 1])))

(define tblis-scalar-set!
  (lambda (ptr x)
    (cond
      [(real? x)
        (ftype-set! tblis-scalar (type) ptr (symbol->tblis-type 'double-float))
        (ftype-set! tblis-scalar (scalar double) ptr (inexact x))])))

(define tblis-scalar->number
  (lambda (ptr)
    (let ([type (ftype-ref tblis-scalar (type) ptr)])
      (cond
        [(equal? type (symbol->tblis-type 'single-float)) (ftype-ref tblis-scalar (scalar float) ptr)]
        [(equal? type (symbol->tblis-type 'double-float)) (ftype-ref tblis-scalar (scalar double) ptr)]))))

(define-syntax foreign-alloc-array
  (syntax-rules ()
    [(_ dtype n) (make-ftype-pointer dtype (foreign-alloc (* (ftype-sizeof dtype) n)))]))

(define-syntax foreign-alloc-object
  (syntax-rules ()
    [(_ dtype) (foreign-alloc-array dtype 1)]))

(define tensor->tblis-tensor
  (lambda (t scale)
    (let* ([ndim (tensor-ndim t)]
           [shape (foreign-alloc-array tblis-len_type ndim)]
           [strides (foreign-alloc-array tblis-stride_type ndim)]
           [ptr (foreign-alloc-object tblis-tensor)])
      (ftype-set! tblis-tensor (type) ptr (symbol->tblis-type (tensor-dtype t)))
      (ftype-set! tblis-tensor (conj) ptr 0)
      (tblis-scalar-set! (ftype-&ref tblis-tensor (scalar) ptr) scale)
      (ftype-set! tblis-tensor (data) ptr (tensor-data t))
      (ftype-set! tblis-tensor (ndim) ptr ndim)
      (do ([i 0 (+ i 1)])
          (= i ndim)
        (ftype-set! tblis-len_type () shape i (tensor-size t i))
        (ftype-set! tblis-stride_type () strides i (tensor-stride t i)))
      (ftype-set! tblis-tensor (len) ptr shape)
      (ftype-set! tblis-tensor (stride) ptr strides)
      ptr)))

(define tblis-tensor-free
  (lambda (ptr)
    (foreign-free (ftype-pointer-address (ftype-ref tblis-tensor (len) ptr)))
    (foreign-free (ftype-pointer-address (ftype-ref tblis-tensor (stride) ptr)))
    (foreign-free (ftype-pointer-address ptr))))

(define with-scaled-tblis-tensor
  (lambda (t scale action)
    (with (delay (tensor->tblis-tensor t scale)) tblis-tensor-free action)))
; (define with-halide-double
;   (lambda (action)
;     (alloca
;       (ftype-sizeof halide_type_t)
;       (lambda (raw)
;         (let ([p (make-ftype-pointer halide_type_t raw)])
;           (ftype-set! halide_type_t (code) p 2)
;           (ftype-set! halide_type_t (bits) p 64)
;           (ftype-set! halide_type_t (lanes) p 1)
;           (action p))))))
(define t (∂-matrix-bounded 0.0 1.0 4))
(display (tensor-get t 0))
(newline)
(tensor-free t)

(do ([i 3 (fx- i 1)])
    ((fx= i 0))
  (malloc-Storage 'double-float 10)
  (collect))
(free-dropped-Storage)
; (display (tensor-get t 1))
; (newline)
; (display (tensor-get t 2))
; (newline)

; (with-halide-double print-halide-type)
; (alloca 10 (lambda (p) (display p) (newline)))
; (define halide_double (halide_type_t halide_type_float 64 1))

