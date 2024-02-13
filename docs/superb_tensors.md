# Quick tutorial

## Tensor creation

There's only one class to declare multidimensional arrays, `Tensor<N,T>`, which the template parameters:

* `N`: number of dimensions, eg `1`, `10` (less than 16, please)
* `T`: datatype of the values, eg `double`. Recommendation: use one of `Complex`, `ComplexF`, `ComplexD`
  for complex values

The `Tensor` class requires labeling all dimensions and give explicit dimensions. For instance,

* Spin matrix $t_{ij}$, Ns x Ns matrix with labels "i" for the first dimension and "j" for the second dimension:
```c++
Tensor<2,Complex> t("ij", {{Ns,Ns}});
```

* Spin-color matrix $t_{cs}$, Nc x Ns matrix with labels "c" for the first dimension and "s" for the second dimension
```c++
Tensor<2,Complex> t("cs", {{Nc,Ns}});
```

* Lattice real field $t_{xyztX}$. Note that `Tensor` deals with the QDP even-odd layout of lattice fields as an extra
  dimension, `X`:
```
Tensor<5,double> t("xyztX", {{Layout::lattSize()[0]/2, Layout::lattSize()[1], Layout::lattSize()[2],
                              Layout::lattSize()[3]}});
// The function `latticeSize` returns the lattice dimensions in the proper order
Tensor<5,double> t("xyztX", latticeSize<5>("xyztX"));
```

* Lattice spin-color field $t_{csxyztX}$:
```
Tensor<7,double> t("csxyztX", latticeSize<7>("csxyztX"));
```

Note that the ordering of the labels indicates how the tensor values arranged in memory, where the shortest stride is for the first dimension.

## Tensor attributes and manipulation

One of the biggest differences with Chroma's objects is that `Tensor` objects are implicit *references* to data memory allocations. For instance,
there's only a single tensor allocation in the following code:
```
Tensor<2,Complex> t("ij", {{Ns,Ns}});
Tensor<2,Complex> q = t; // t and q refer the same allocation
```
The equivalent Chroma code is:
```
SpinMatrix *t = new SpinMatrix;
SpinMatrix *q = t;
```
The memory allocations are automatically freed when they don't have references.
(If you know what a smart pointer is, then `Tensor` has a `std::shared_ptr` to an allocation.)

List of basic tensor attributes:

* `order`, the dimension labels.
* `dim`, the tensor dimensions.
* `from`, first element coordinates in the allocation.
* `size`, dimensions of these tensor.
* `volume()`, return the number of elements
* `operator bool()`, return whether the volume isn't zero.
* `kvdim() -> std::map<char,int>`, return the size associated to each dimension label.
* `get(Coor<N> coor) -> T`, return the value of tensor at that coordinate.

List of basic tensor operations:

* `rename_dims(std::map<char,char> remap) -> Tensor<N,T>`, return a view with the dimension labels renamed.
* `kvslice_from_size(std::map<char,int> from, std::map<char,int> size) -> Tensor<N,T>`, returns a slice of the current view.
* `like_this<Nn=N,Tn=T>(Maybe<std::string> order, std::map<char,int> size={}) -> Tensor<Nn,Tn>`,
  return a new allocation
* `clone() -> Tensor<T,N>`, return a new allocation with the same values.
* `reorder(std::string order) -> Tensor<T,N>`, return a new allocation with the elements ordered as given; it may return the same allocation.
* `set_zero()`, set to zero all tensor elements.
* `copyTo(Tensor<Nn,Tn> t)`, copy all tensor elements into the given tensor.
* `release()`, remove the reference to the allocation.

The following code shows the correspondence between Chroma objects and `Tensor` objects:
```C++
// Chroma objects                   | // New objects
SpinMatrix t;                       | Tensor<2,Complex> t("ij", {{Ns,Ns}});
                                    | t.set_zero(); // Tensor objects are uninitialized!
SpinMatrix q = t;                   | Tensor<2,Complex> q = t.clone();

LatticeSpin x;                      | Tensor<6,Complex> x("sxyzXt", latticeSize("sxyzXt"));
LatticeSpin y = x * Gamma(1);       | Tensor<6,Complex> y = x.like_this();
                                    | y.contract(x, {{'x','i'}}, NotConjugate,
                                    |            asTensorView(Gamma(i)), {}, NotConjugate, {{'j','s'}});
```

## Tensor arithmetic operations

* Scalar
```
Tensor<2,Complex> t("ij", {{Ns,Ns}});
auto q = t.scale(-1); // q_ij = -t_ij
```

* Addition:
```
Tensor<2,Complex> t("ij", {{Ns,Ns}});
Tensor<2,Complex> q("ij", {{Ns,Ns}});
t.addTo(q); // q_ij += t_ij
```

* Contraction, contract the labels that are on the input tensors and not on the output tensor:
```
Tensor<2,Complex> t("ij", {{Ns,Ns}});
Tensor<2,Complex> q("ik", {{Ns,Nc}});
Tensor<2,Complex> r("ik", {{Ns,Nc}});
r.contract(t, {}, NotConjugate, q, {}, NotConjugate); // r_ik = \sum_i t_ij * q_jk
```


