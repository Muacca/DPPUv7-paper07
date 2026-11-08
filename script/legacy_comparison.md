# Conventions and the previous local block

Sections 2 and Appendices A–B specify the current action, fixed quotients and physical radii. The former torsion input array is $`\tau=-T`$ for geometric torsion, and the NY notation converts as $`\theta=\theta_{\mathrm{NY}}/\kappa^2`$. A spherical structure magnitude a/r gives q=2r/|a|, including a=1,2,4. Radius, volume and wave-number units must be converted together.

The trace block (D-04) in Appendix I uses diagonal coframe perturbations s/3, lapse perturbation ℓ, connection boost diagonal entries b, rotation entries 12_3 and 23_1 equal to a, entry 13_2 equal to −a, and other entries zero. Matching the boundary term replaces −8s dot b by 8 dot s b. Eliminating a=−s/6−ℓ/2 and b=−dot s/3, then ℓ=s/3, gives L=−4 dot s²/3.

This strictly quadratic block has 0 first-class constraints, 6 second-class constraints and 1 pair in the connection description, or 0 first-class constraints, 2 second-class constraints and 1 pair after elimination. The original 26-configuration block has 0 first-class constraints, 50 second-class constraints and 1 pair. Primary preservation includes the multiplier term of the primary Poisson matrix, rather than setting every primary–Hamiltonian bracket to zero. Perturbation lapse times connection squared is separated as a cubic term in the previous local action.

This block is distinct from the scalar quadratic action derived from the complete action in Appendix C. In `routes.json`, `legacy_trace` is this comparison; `0` and `1` are the complete-action torus and sphere controls. Appendices D and F–H use the same action conventions. See the [changelog](../CHANGELOG.md) for changes between editions and the [script guide](README.md) for current reproduction instructions.
