#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Placeholder model.

FedAvgRecipe requires a model source, but this job does not train anything:
the summary statistics travel in FLModel.meta and the meta-analysis happens in
the aggregator. This exists only to satisfy the recipe and to give the model
persistor a real state dict to write.
"""

import torch
import torch.nn as nn


class DummyModel(nn.Module):
    """A single parameter, carried round-trip and otherwise unused."""

    def __init__(self):
        super().__init__()
        self.placeholder = nn.Parameter(torch.zeros(1))

    def forward(self, x):
        return x