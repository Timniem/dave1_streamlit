import plotly.graph_objects as go

def explain_plot_plotly(shap_values, raw_values, feature_names, predicted_score):
    # Prepare data
    contributions = [(name, val) for name, val in zip(feature_names, shap_values)]
    contributions = sorted(contributions, key=lambda x: abs(x[1]))

    cumulative = predicted_score
    bars = []
    bases = []
    hover_texts = []
    colors = []
    y_labels = []

    for name, val in contributions:
        neg_val = -val
        bars.append(neg_val)
        bases.append(cumulative)
        hover_texts.append(f"{val:.3f}")
        colors.append("#FF0C57" if val > 0 else "#017FFD")
        if raw_values[name] % 1 == 0:
            y_labels.append(f"{name} ({raw_values[name]})")
        else:
            y_labels.append(f"{name} ({raw_values[name]:.03f})")
        cumulative += neg_val

    fig = go.Figure()

    fig.add_trace(go.Bar(
        x=bars,
        y=y_labels,
        orientation='h',
        base=bases,
        marker=dict(color=colors),
        hovertext=hover_texts,
        hoverinfo="text",
        showlegend=False
    ))

    fig.update_layout(
        xaxis=dict(
        title="DAVE1 score",
        mirror='ticks',
        showgrid=True,
        range=[0, 1],
        ),
        yaxis=dict(
            showticklabels=True,
            showgrid=True,
        ),
        plot_bgcolor='white',
        margin=dict(l=20, r=20, t=0, b=20),
    )

    fig.layout.xaxis.fixedrange = True
    fig.layout.yaxis.fixedrange = True
    

    return fig
