import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np


def plot_scores(shap_scores, raw_scores, ax, predicted_score, maximum_value):
    cumulative = [predicted_score,]
    for n,i in enumerate(shap_scores):
        score = shap_scores[i]
        color = "#FF0C57" if score > 0 else "#1E88E5"
        
        neg_score = -score
        if n + 1 < len(shap_scores):
            cumulative.append(cumulative[n] + neg_score)

        bar_start = cumulative[n] + neg_score / 2
        ax.barh(0, neg_score, left=cumulative[n], height=0.5, color=color, edgecolor='white', alpha=1, linewidth=1.5)
        points = [(bar_start, 0), (cumulative[n] + neg_score, -0.25), (cumulative[n] + neg_score, 0.25)]
        polygon = Polygon(points, closed=True, edgecolor='white', facecolor='white', alpha=0.4, linewidth=0)

        ax.add_patch(polygon)

        if abs(score) < maximum_value / 10 or n >= 2:
            continue
            
        if n % 2 == 0:
            orientation = 1
        else:
            orientation = -1
        # Draw vertical line from bar to annotation
        ax.plot([bar_start, bar_start], [-.3* orientation, -0.5 * orientation], color='#BBB', linestyle='-', linewidth=2)
        # Place annotation below the plot
        ax.text(bar_start, -0.75 * orientation, f'{i} ({raw_scores[i]:.3f}) = {score:.3f}', ha='center', va='top', fontsize=14, color='#444')


def force_plot(shap_values, raw_values, feature_names, predicted_score):

    # Separate positive and negative values
    positives = dict(sorted({nom:num for num,nom in zip(shap_values, feature_names) if num > 0}.items(), key=lambda item: item[1], reverse=True))
    negatives = dict(sorted({nom:num for num,nom in zip(shap_values, feature_names) if num < 0}.items(), key=lambda item: item[1]))

    if len(negatives) > 0 and len(positives) > 0:
        maximum_value = max(positives[next(iter(positives))], abs(negatives[next(iter(negatives))]))
    elif len(negatives) == 0 and len(positives) > 0:
        maximum_value = positives[next(iter(positives))]
    else:
        maximum_value = abs(negatives[next(iter(negatives))])

    # Plot setup
    fig, ax = plt.subplots(figsize=(12, 3))  # Increased height for annotations

    plot_scores(positives, raw_values, ax, predicted_score, maximum_value)
    plot_scores(negatives, raw_values, ax, predicted_score, maximum_value)

    ax.plot([predicted_score, predicted_score], [-.3, .3], color="#F1F1F1", linestyle='-', linewidth=5)
    ax.text(predicted_score, -0.65, f'{predicted_score:.3f}', ha='center', va='bottom', fontsize=20, color='#000', fontweight='bold')
    
    # Formatting
    ax.set_yticks([])
    ax.xaxis.set_ticks_position('top') # the rest is the same

    ax.set_ylim(-1, 1)  # Adjusted to accommodate annotations
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(color='#EEE', labelcolor='#333')
    ax.spines['top'].set_edgecolor('#EEE')
    fig.patch.set_facecolor('#FFFFFF')
    ax.set_facecolor('#FFFFFF')
    plt.tight_layout()
    return fig