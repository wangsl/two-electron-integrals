

      Program Main
      Implicit None
      Integer N, I
      Parameter(N=5)
      Real*8 r(N), w(N)
      Real*8 x

      Call RysChebyshevInitialize
      
      x = 0.0D0
      Do while(x <= 10.000001D0)
        Call CalculateRysRootsAndWeights(n, x, r, w, 0)
        Write(6, *)
        Do I = 1, n
          Write(6, '( F8.2  I4 F22.18 E26.18)') x, I, r(I), w(I)
        EndDo
        x = x + 1.0D0
      EndDo

      End