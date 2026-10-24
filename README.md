# Problem 17 — PIN search with Grover's algorithm

Finds the hidden 4-bit PIN `1011` among 16 possibilities using Qiskit and the Aer simulator.

```
pip install -r requirements.txt
python grover_pin.py
```

Outputs: `grover_histogram.png` (3 iterations, 1024 shots) and `grover_success_vs_iterations.png` (simulated vs theory, k = 1..6).

## Explanation

1. **Oracle:** X gates on the qubits where the PIN has a `0` map |1011⟩ to |1111⟩. A multi-controlled Z (H·MCX·H) flips its phase, then the X gates are undone, so only |1011⟩ picks up a −1. Qiskit puts qubit 0 at the rightmost position of a bitstring, so the PIN is reversed before choosing qubits. That makes the measured `1011` really mean 1011.
2. **Diffuser:** H → X → multi-controlled Z → X → H applies 2|s⟩⟨s| − I, which reflects every amplitude about the mean. This amplifies the marked state that the oracle made negative.
3. **Theory:** with θ = arcsin(1/√16) = arcsin(1/4), after k iterations P(1011) = sin²((2k+1)θ). The optimum is k = ⌊π/4·√16⌋ = 3, giving P ≈ 0.961.
4. **Results:** with 3 iterations, `1011` was measured 997/1024 times (97.4%). Simulated vs theory for k = 1..6: 0.502/0.473, 0.929/0.908, 0.974/0.961, 0.604/0.582, 0.121/0.126, 0.020/0.020. Every value is within 0.03 of theory, which is shot noise for 1024 shots.
5. **Overshoot:** after k = 3 the state rotates past the target, so the probability drops (≈0.58 at k = 4, ≈0.13 at k = 5). Too many iterations hurts, which is why Grover needs the iteration count chosen in advance.
