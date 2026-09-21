"""Portable image-comparison helpers for plot regression tests."""

import numpy as np
from PIL import Image, ImageFilter
from skimage.metrics import structural_similarity


# Font rendering and tight bounding boxes vary across OS/Matplotlib versions.
# Compare at a common resolution, but reject large canvas-size changes before
# resizing so a genuinely different layout cannot be normalized away.
MINIMUM_STRUCTURAL_SIMILARITY = 0.94
MINIMUM_CHROMA_SIMILARITY = 0.95
IMAGE_BLUR_RADIUS = 1
MINIMUM_SATURATION = 0.08
MAXIMUM_DIMENSION_VARIANCE = 0.10
MAXIMUM_COMPARISON_WIDTH = 1400


def _blurred_rgb(image):
    return image.convert("RGB").filter(ImageFilter.GaussianBlur(IMAGE_BLUR_RADIUS))


def _chroma_descriptor(image):
    """Describe foreground color independently of white-background area."""
    hsv = np.asarray(image.convert("HSV"), dtype=np.float64)
    hue = hsv[:, :, 0]
    saturation = hsv[:, :, 1] / 255.0
    value = hsv[:, :, 2] / 255.0
    foreground = (saturation >= MINIMUM_SATURATION) & (value < 0.995)
    weights = saturation * foreground
    chroma_mass = float(weights.sum() / weights.size)
    histogram = np.histogram(
        hue, bins=24, range=(0, 256), weights=weights
    )[0]
    if histogram.sum():
        histogram = histogram / histogram.sum()
    return chroma_mass, histogram


def _chroma_similarity(candidate, reference):
    candidate_mass, candidate_histogram = _chroma_descriptor(candidate)
    reference_mass, reference_histogram = _chroma_descriptor(reference)

    if not candidate_mass and not reference_mass:
        return 1.0
    if not candidate_mass or not reference_mass:
        return 0.0

    mass_similarity = min(candidate_mass, reference_mass) / max(
        candidate_mass, reference_mass
    )
    histogram_intersection = float(
        np.minimum(candidate_histogram, reference_histogram).sum()
    )
    return min(mass_similarity, histogram_intersection)


def image_similarities(candidate_path, reference_path):
    """Return structural and foreground-chroma similarity scores."""
    with Image.open(candidate_path) as candidate, Image.open(reference_path) as reference:
        if any(
            abs(candidate_extent - reference_extent) / reference_extent
            > MAXIMUM_DIMENSION_VARIANCE
            for candidate_extent, reference_extent in zip(candidate.size, reference.size)
        ):
            return 0.0, 0.0
        candidate = candidate.convert("RGB").resize(
            reference.size, Image.Resampling.LANCZOS
        )
        if reference.width > MAXIMUM_COMPARISON_WIDTH:
            comparison_size = (
                MAXIMUM_COMPARISON_WIDTH,
                round(reference.height * MAXIMUM_COMPARISON_WIDTH / reference.width),
            )
            candidate = candidate.resize(comparison_size, Image.Resampling.LANCZOS)
            reference = reference.resize(comparison_size, Image.Resampling.LANCZOS)
        candidate = _blurred_rgb(candidate)
        reference = _blurred_rgb(reference)
        structure = structural_similarity(
            np.asarray(candidate),
            np.asarray(reference),
            data_range=255,
            channel_axis=2,
        )
        return structure, _chroma_similarity(candidate, reference)
