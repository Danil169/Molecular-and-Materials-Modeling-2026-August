Anharmocity of water with MD
============================


https://share.google/aimode/sm3TizyMqtQMZ9o4d


working code
~~~~~~~~~~~~
(venv) milias@DESKTOP-7OTLCGO:~/work/projects/Molecular-and-Materials-Modeling-2026-May/water-MD-anharmonicity/.python code4.py
CHGNet v0.3.0 initialized with 412,525 parameters
CHGNet will run on cuda
--- Equilibration at 300 K ---
/home/milias/work/software/venv/lib/python3.12/site-packages/chgnet/model/model.py:898: UserWarning: Converting a tensor with requires_grad=True to a scalar may lead to unexpected behavior.
Consider using tensor.detach() first. (Triggered internally at /pytorch/torch/csrc/autograd/generated/python_variable_methods.cpp:836.)
  volumes = torch.tensor(volumes, dtype=TORCH_DTYPE, device=atomic_numbers.device)
--- Production Run ---



Results at Temperature: 300 K
--------------------------------------------------
Bending        : Calc  1622.6 | Expt  1595.0 | Error  1.73%
Sym-Stretch    : Calc  3652.4 | Expt  3657.0 | Error  0.12%
Asym-Stretch   : Calc  3754.2 | Expt  3756.0 | Error  0.05%
