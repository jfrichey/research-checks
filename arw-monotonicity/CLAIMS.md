# Proof and verification map

## Proof records

| Reader key | Supporting argument |
| --- | --- |
| RN | [Response renewal](response-renewal.md), §§1–37 |
| DI | [Drain inventories](drain-inventory.md), §§1–12 |
| SX | [Stack exploration](stack-exploration.md) |
| CX | [Creator exchanges](creator-exchange.md) |
| AX | [Avalanche exchanges](avalanche-exchange.md) |
| NA-report | [Sink reinjection](sink-reinjection.md), (50)–(58) |
| Earlier letter keys | [Evidence map](evidence/README.md) |

The reader supplies the section/equation with each source key. These linked arguments retain their working, incomplete or reported status; an AI audit is not independent human review.

## Exact verification map

| Family in the reader | Program / saved certificate |
| --- | --- |
| Local word and path exchanges | [`verify_stack_exploration.py`](verify_stack_exploration.py); [`verify_drained_paths.py`](verify_drained_paths.py) |
| Response and creation laws | [`verify_response_renewal.py`](verify_response_renewal.py); [`verify_creation_cone.py`](verify_creation_cone.py); [`verify_root_restoration.py`](verify_root_restoration.py) |
| Creator potentials and contacts | [`verify_creator_potentials.py`](verify_creator_potentials.py); [`verify_contact_renewal.py`](verify_contact_renewal.py); [`verify_backward_kernels.py`](verify_backward_kernels.py) |
| Three-site occupation | [`verify_three_site_certificate.py`](verify_three_site_certificate.py); [`three-site-positive-certificate.json.gz`](three-site-positive-certificate.json.gz); [`verify_three_site_renewal.py`](verify_three_site_renewal.py) |
| Pair deficit capacity | [`verify_three_site_pair_deficit.py`](verify_three_site_pair_deficit.py); [`three-site-pair-deficit-certificate.json.gz`](three-site-pair-deficit-certificate.json.gz); [`verify_pair_capacity.py`](verify_pair_capacity.py) |
| Inventory positivity | [`verify_drain_inventory.py`](verify_drain_inventory.py); [`verify_single_drain_inventory.py`](verify_single_drain_inventory.py); [`verify_inventory_corners.py`](verify_inventory_corners.py); [`verify_cycle_inventory.py`](verify_cycle_inventory.py); [`verify_three_site_permanent.py`](verify_three_site_permanent.py) |
| Directed queues and feedback | [`verify_cycle_queue.py`](verify_cycle_queue.py); [`verify_cycle_renewal.py`](verify_cycle_renewal.py); [`verify_backward_queue_cone.py`](verify_backward_queue_cone.py); [`verify_feedback_cone.py`](verify_feedback_cone.py) |
| Small-sleep and boundary calculations | [`verify_three_site_rare_sleep.py`](verify_three_site_rare_sleep.py); [`verify_cycle_rare_sleep.py`](verify_cycle_rare_sleep.py); [`verify_three_site_boundary.py`](verify_three_site_boundary.py) |
| Sleep transport and resident release | [`verify_sleep_transport.py`](verify_sleep_transport.py); [`verify_release_forest.py`](verify_release_forest.py); [`verify_cycle_return.py`](verify_cycle_return.py) |
| Sink reinjection and stationary rank | [`weak_interaction.py`](../arw-negative-association/weak_interaction.py); [`exact.py`](../arw-negative-association/exact.py) |
| Markers and filling structures | [`verify_marker_geometry.py`](verify_marker_geometry.py); [`verify_sleep_valuations.py`](verify_sleep_valuations.py) |
| Boundary responses and excursions | [`verify_boundary_response.py`](verify_boundary_response.py); [`verify_excursion_dual.py`](verify_excursion_dual.py) |
| Currents and addition geometry | [`verify_exit_tilt.py`](verify_exit_tilt.py); [`verify_addition_algebra.py`](verify_addition_algebra.py) |

The default suite runs `verify_cycle_return.py`: six symbolic identities and 2,232 Bernstein coefficients for Theorem 8.1. Other programs preserve the previously reported exact audits and require their stated modes and resources. They were not all rerun for this editorial release. The response-renewal loader now reads two selected archive members directly; its default audit was rerun to check that packaging change. Saved compressed certificates are identified in `../files.sha256.json`; checking a file hash is distinct from regenerating its symbolic proof.

Use Python 3.12 with the root `requirements-checks.txt`; optional large symbolic modes may also require `python-flint`. Run programs from this directory to preserve local input paths. Never use Python `-O`, which disables their mathematical assertions.
