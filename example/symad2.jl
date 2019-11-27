using Yao, YaoExtensions
@vars α β γ
c = chain(put(3, 2=>Rx(α)), control(3, 2, 1=>Ry(β)), put(3, (1,2)=>rot(kron(X, X), γ)))
h = heisenberg(3)
expect(h, zero_state(Basic, 3)=>c)
expect'(h, zero_state(Basic, 3)=>c)
