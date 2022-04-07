import plotly.graph_objects as go
from plotly.offline import iplot
import pandas as pd
import numpy as np
import plotly.express as px


def plot_single_clf(graph_metrics, quality_metrics, title):
    fig = go.Figure()
    prox_metric = list(graph_metrics.keys())
    clf_dfs = list(graph_metrics.values())

    for metric in quality_metrics:
        xlabels = [metric.capitalize()]*len(clf_dfs[0].index)*len(clf_dfs)
        xlabels2 = [prox_metric[0]]*int(len(xlabels)/len(clf_dfs))
        y_data = clf_dfs[0][metric].values
        for i in range(1, len(clf_dfs)):
            xlabels2 += [prox_metric[i]]*int(len(xlabels)/len(clf_dfs))
            y_data = np.concatenate(
                [y_data, clf_dfs[i][metric].values], axis=0)
        xlabels1 = [xlabels2, xlabels]
        fig = fig.add_trace(go.Box(y=y_data, showlegend=False,
                            boxmean=True, width=0.6, x=xlabels1, jitter=0.8))
    fig.update_layout(
        title=title,
        yaxis_title="Score",
        xaxis_title="Metric/Network",
        boxmode='group',
        yaxis=dict(autorange=True, showgrid=True, zeroline=True, dtick=0.1, gridcolor='rgb(255, 255, 255)',
                   gridwidth=1, zerolinecolor='rgb(255, 255, 255)', zerolinewidth=2),
        title_x=0.5, paper_bgcolor='rgb(255, 255, 255)', plot_bgcolor='rgb(243, 243, 243)')

    fig.show()


def plot_clf_comparision(graph_metrics1, graph_metrics2, clf_type, quality_metric, title):
    fig = go.Figure()
    prox_metric = list(graph_metrics1.keys())
    clf_dfs1 = list(graph_metrics1.values())
    clf_dfs2 = list(graph_metrics2.values())

    for clf in clf_type:
        if clf == clf_type[0]:
            xlabels = [clf]*len(clf_dfs1[0].index)*len(clf_dfs1)
            xlabels2 = [prox_metric[0]]*int(len(xlabels)/len(clf_dfs1))
            y_data = clf_dfs1[0][quality_metric].values
            for i in range(1, len(clf_dfs1)):
                xlabels2 += [prox_metric[i]]*int(len(xlabels)/len(clf_dfs1))
                y_data = np.concatenate(
                    [y_data, clf_dfs1[i][quality_metric].values], axis=0)
            xlabels1 = [xlabels2, xlabels]
        else:
            xlabels = [clf]*len(clf_dfs2[0].index)*len(clf_dfs2)
            xlabels2 = [prox_metric[0]]*int(len(xlabels)/len(clf_dfs2))
            y_data = clf_dfs2[0][quality_metric].values
            for i in range(1, len(clf_dfs1)):
                xlabels2 += [prox_metric[i]]*int(len(xlabels)/len(clf_dfs1))
                y_data = np.concatenate(
                    [y_data, clf_dfs2[i][quality_metric].values], axis=0)
            xlabels1 = [xlabels2, xlabels]
        fig = fig.add_trace(go.Box(y=y_data, showlegend=False,
                            boxmean=True, width=0.6, x=xlabels1, jitter=0.8))

    fig.update_layout(
        title=title,
        yaxis_title=quality_metric.capitalize(),
        xaxis_title="Metric/Network",
        boxmode='group',
        yaxis=dict(
            autorange=True,
            showgrid=True,
            zeroline=True,
            dtick=0.1,
            gridcolor='rgb(255, 255, 255)',
            gridwidth=1,
            zerolinecolor='rgb(255, 255, 255)',
            zerolinewidth=2), title_x=0.5, paper_bgcolor='rgb(255, 255, 255)', plot_bgcolor='rgb(243, 243, 243)')

    fig.show()


def multiview(prox_metrics_list, metrics):

    fig = go.Figure()
    buttons = []
    i = 0
    colors = px.colors.qualitative.Plotly[:3]
    # add line and shaded area for each series and standards deviation

    for prox_metric_label, prox_metric in prox_metrics_list.items():
        for label_, data in prox_metric.items():

            for pos, metric in enumerate(metrics):
                xlabels = [label_.capitalize()]*len(data)
                xlabels2 = [metric]*int(len(xlabels)/len(data))
                y_data = data[metric].values
                for z in range(1, len(prox_metric)):
                    #xlabels2 += [metric]*int(len(xlabels)/len(data))
                    xlabels2 += [metric]*int(len(xlabels)/len(metrics))
                    y_data = np.concatenate(
                        [y_data, data[metric].values], axis=0)
                    #print(xlabels2)
                xlabels1 = [xlabels, xlabels2]
                # print(xlabels1)
                fig = fig.add_trace(go.Box(y=y_data, showlegend=False, boxmean=True,
                                    width=0.6, x=xlabels1, jitter=0.8, marker_color=colors[pos]))

            # args is a list of booleans that tells the buttons which trace to show on click
        args = [False] * len(prox_metric)*len(metrics)*len(prox_metrics_list)
        n_box = len(metrics)*len(prox_metric)
        args[i*n_box:(i+1)*n_box] = [True]*n_box
        # create a button object for the country we are on
        button = dict(label=prox_metric_label.replace('_', ' ').title(),
                      method="update",
                      args=[{"visible": args}])
        # add the button to our list of buttons
        buttons.append(button)
        # i is an iterable used to tell our "args" list which value to set to True
        i += 1

    fig.update_layout(
        updatemenus=[
            dict(
                # change this to "buttons" for individual buttons
                type="buttons",
                # this can be "left" or "right" as you like
                # direction="down",
                # (1,1) refers to the top right corner of the plot
                x=0,
                xanchor='left',
                y=1.01,
                yanchor='bottom',
                direction='left',
                # the list of buttons we created earlier
                buttons=buttons)
        ])
    # fig.update_layout(
    # title='title',
    # yaxis_title="Score",
    # xaxis_title="Metric/Network",
    # boxmode='group',
    # yaxis=dict(autorange=True, showgrid=True, zeroline=True, dtick=0.1, gridcolor='rgb(255, 255, 255)',
    #           gridwidth=1, zerolinecolor='rgb(255, 255, 255)', zerolinewidth=2),
    #           title_x=0.5, paper_bgcolor='rgb(255, 255, 255)', plot_bgcolor='rgb(243, 243, 243)')
    fig.update_layout(yaxis=dict(range=[-0.2, 1.2], ), title='Classification Results Across Different Metrics',
                      xaxis_title='Process Number', yaxis_title='F-Measure', title_xanchor='left', title_x=0.4)
    fig.show()
    return
