# print.kernel_spec is stable

    Code
      print(kernel_spec())
    Output
      Kernel specification:
        Type: rbf 
        Bandwidth: median heuristic 

---

    Code
      print(kernel_spec("rbf", bandwidth = 1.5))
    Output
      Kernel specification:
        Type: rbf 
        Bandwidth: 1.5 

---

    Code
      print(kernel_spec("matern", bandwidth = 1, nu = 1.5))
    Output
      Kernel specification:
        Type: matern 
        Bandwidth: 1 
        nu: 1.5 

---

    Code
      print(kernel_spec("rbf", approx = "nystrom", approx_rank = 50, approx_seed = 1))
    Output
      Kernel specification:
        Type: rbf 
        Bandwidth: median heuristic 
        Approx: nystrom 
        Approx rank: 50 
        Approx seed: 1 

---

    Code
      print(kernel_spec("rbf", approx = "rff", approx_rank = 200, approx_seed = 7))
    Output
      Kernel specification:
        Type: rbf 
        Bandwidth: median heuristic 
        Approx: rff 
        Approx rank: 200 
        Approx seed: 7 

# print.nystrom_factor is stable

    Code
      print(nf)
    Output
      Nystrom factor:
        n        : 20 
        m        : 10 
        rank     : 10 
        kernel   : rbf 

# print.rff_features is stable

    Code
      print(rf)
    Output
      Random Fourier features:
        n        : 20 
        D        : 25 
        features : 50 (2 * D)
        kernel   : rbf 

