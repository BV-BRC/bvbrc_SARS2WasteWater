#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 06:33:08 2024

@author: nbowers
"""
import ast  # Import ast module for literal_eval
import matplotlib.colors as mcolors
import pandas as pd
import plotly.graph_objects as go
import sys


color_blind_friendly_colors = [
    '#E69F00',  # Orange
    '#56B4E9',  # Sky Blue
    '#009E73',  # Bluish Green
    '#F0E442',  # Yellow
    '#0072B2',  # Blue
    '#D55E00',  # Vermilion
    '#CC79A7',  # Reddish Purple
    '#999999',  # Gray
]

def get_extended_color_palette(n_colors):
    return [mcolors.hsv_to_rgb((x*1.0/n_colors, 0.5, 0.9)) for x in range(n_colors)]

def get_lineage_info(df):
    lineage_df = df[["sample", "lineages", "abundances"]]

    # Split the 'lineages' and 'abundances' columns into lists
    lineage_df['lineages'] = df['lineages'].str.split()
    lineage_df['abundances'] = df['abundances'].str.split().apply(lambda x: [float(i) for i in x])

    # Expand the lists into rows
    rows = []
    for _, row in lineage_df.iterrows():
        for lineage, abundance in zip(row['lineages'], row['abundances']):
            rows.append({'sample': row['sample'], 'lineage': lineage, 'abundance': abundance})

    # Create a new DataFrame
    df_lineages = pd.DataFrame(rows)
    # drop lineage from plotting if it is less than 2%
    df_lineages = df_lineages[df_lineages['abundance'] >= 0.02]
    return df_lineages

def get_variant_info(df):
    # Subsection only the data we need
    variants_df = df[["sample", "summarized"]]

    variants_df['summarized'] = variants_df['summarized'].apply(ast.literal_eval)
    # Parse the 'Variants' column to a more manageable form
    rows = []
    for index, row in variants_df.iterrows():
        for variant, percent in row['summarized']:
            rows.append({'Sample': row['sample'], 'Variant': variant, 'Percent': percent})
    df_variants = pd.DataFrame(rows)
    return df_variants

def plot_lineage_by_samples(df_lineages, sample_lineage_out):
    # Create a Plotly figure
    fig = go.Figure()
    # # changes to atch below 
    # Add traces
    lineages = df_lineages['lineage'].unique()
    for i, lineage in enumerate(lineages):
        if len(lineage) <= len(color_blind_friendly_colors):
            color = color_blind_friendly_colors[i % len(color_blind_friendly_colors)]
        else:
            n = len(lineage)
            color = get_extended_color_palette(n)
        df_lineage = df_lineages[df_lineages['lineage'] == lineage]
        fig.add_trace(go.Bar(
            x=df_lineage['sample'],
            y=df_lineage['abundance'],
            name=lineage,
            marker_color=color,
            hoverinfo='y+name',
            hovertemplate='<b>%{x}</b><br>%{y:.2f}%<br><b>%{data.name}</b><extra></extra>'  # Customize hovertext with HTML
            # hoverlabel=dict(font_color='black')  # Set hover text color to black
        ))

    # Update layout for a stacked bar chart
    fig.update_layout(
        barmode='stack',
        title='Lineage Abundances per Sample',
        title_font_size=24,
        xaxis_title='Sample',
        xaxis_title_font_size=18,
        yaxis_title='Abundance (%)',
        yaxis_title_font_size=18,
        yaxis=dict(tickformat=".0%", tickfont_size=16),
        legend_title='Lineage',
        legend_title_font_size=16,
        legend_font_size=14,
        xaxis=dict(tickangle=-45, tickfont_size=16),
        hoverlabel=dict(font_size=16, font_family="Arial")
    )

    # Display the figure or save it as HTML
    fig.write_html(sample_lineage_out)  # Save the interactive plot as an HTML file
    return


def plot_variant_by_samples(df_variants, sample_variant_out):
    # Create a Plotly figure
    fig = go.Figure()
    # # changes to atch below 
    # Add traces
    variants = df_variants['Variant'].unique()
    for i, variant in enumerate(variants):
        df_subset = df_variants[df_variants['Variant'] == variant]
        color = color_blind_friendly_colors[i % len(color_blind_friendly_colors)]
        fig.add_trace(go.Bar(
            x=df_subset['Sample'],
            y=df_subset['Percent'],
            name=variant,
            marker_color=color,  # Set color from the palette
            hoverinfo='y+name',
            hovertemplate='<b>%{x}</b><br>%{y:.2f}%<br><b>%{data.name}</b><extra></extra>'  # Customize hovertext with HTML
        ))

    # Update layout for a stacked bar chart
    fig.update_layout(
        barmode='stack',
        title='Variant Percentages per Sample',
        title_font_size=24,
        xaxis_title='Sample',
        xaxis_title_font_size=18,
        yaxis_title='Percentage',
        yaxis_title_font_size=18,
        yaxis=dict(tickformat=".0%", tickfont_size=16),
        legend_title='Variant',
        legend_title_font_size=16,
        legend_font_size=14,
        xaxis=dict(tickangle=-45),  # Rotate labels to -45 degrees
        hoverlabel=dict(font_size=16, font_family="Arial")
    )


    # Display the figure in HTML or save it
    fig.write_html(sample_variant_out)  # This saves the interactive plot as an HTML file
    return

def main(argv):
    # step 0 get paths set up
    freyja_results = argv[1]
    sample_lineage_out = argv[2]
    sample_variant_out = argv[3]
    print(sample_lineage_out)
    print(sample_variant_out)

    df = pd.read_csv(freyja_results, sep="\t")
    df.columns = ["sample", "summarized", "lineages", "abundances", "resid","coverage"]
    # chop off the file name and extension 
    df["sample"] = df["sample"].apply(lambda x: x.replace('_freyja_variants.tsv', ''))
    # lineages
    df_lineages = get_lineage_info(df)
    plot_lineage_by_samples(df_lineages, sample_lineage_out)
    # variants
    df_variants = get_variant_info(df)
    plot_variant_by_samples(df_variants, sample_variant_out)
    print("Stacked bar plots generated ")

if __name__ == "__main__":
    main(sys.argv)
