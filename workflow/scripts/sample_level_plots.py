#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 06:33:08 2024

@author: nbowers
"""
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import sys

colorblind_palette = [
    "#377eb8",  # Blue
    "#ff7f00",  # Orange
    "#4daf4a",  # Green
    "#f781bf",  # Pink
    "#a65628",  # Brown
    "#984ea3",  # Purple
    "#999999",  # Grey
    "#e41a1c",  # Red
    "#dede00",  # Yellow
    "#e7298a",  # Dark Pink
    "#66a61e",  # Light Green
    "#1f78b4",  # Light Blue
    "#a6cee3",  # Pale Blue
    "#fdbf6f",  # Light Orange
    "#b2df8a"   # Light Green
]

def get_extended_color_palette(n_colors):
    return [mcolors.hsv_to_rgb((x*1.0/n_colors, 0.5, 0.9)) for x in range(n_colors)]

def get_lineage_info(df):    
    lineage_df = df[["sample", "lineages", "abundances"]].copy()
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
    rows = []
    # Subsection only the data we need
    samples_variants_info = df[["sample", "summarized"]]
    # Will use the entire variant name
    for index, row in samples_variants_info.iterrows():
        sample_id = row["sample"]
        tuples_list = row['summarized']
        tuples_list = list(eval(tuples_list))
        for tup in tuples_list:
            name, abundance = tup
            # Append a new row to the rows list
            rows.append({'sample_id': sample_id, 'Variant': name, 'abundance': abundance})
        # Create a new dataframe from the rows list
    variants_df = pd.DataFrame(rows)
    return variants_df

def plot_lineage_by_samples(df_lineages, sample_lineage_out):
    # Create a Plotly figure
    fig = go.Figure()
    # Add traces
    lineages = df_lineages['lineage'].unique()
    for i, lineage in enumerate(lineages):
        df_lineage = df_lineages[df_lineages['lineage'] == lineage]
        fig.add_trace(go.Bar(
            x=df_lineage['sample'],
            y=df_lineage['abundance'],
            name=lineage,
            marker_color=colorblind_palette[i % len(colorblind_palette)], # loop through the color palette
            hoverinfo='y+name',
            hovertemplate='<b>%{x}</b><br>%{y:.2%}<br><b>%{data.name}</b><extra></extra>',
    ))

    # Update layout for a stacked bar chart
    fig.update_layout(
        barmode='stack',
        title='Lineage Abundance per Sample',
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
        hoverlabel=dict(font_size=16, font_family="Roboto"),
        height=700,
    )
    fig.write_html(sample_lineage_out, include_plotlyjs=False)  # This plot will not work outside of the report
    return


def plot_variant_by_samples(df_variants, sample_variant_out):
    # Create a list of unique variants
    variants = df_variants['Variant'].unique()
    # Initialize a figure
    fig = go.Figure()

    # Add a bar for each variant
    for i, variant in enumerate(variants):
        variant_data = df_variants[df_variants['Variant'] == variant]
        fig.add_trace(go.Bar(
            x=variant_data['sample_id'],
            y=variant_data['abundance'],
            marker_color=colorblind_palette[i % len(colorblind_palette)], # loop through the color palette
            name=variant
        ))

    fig.update_layout(
    barmode='stack',
    title='Variant Abundance by Sample',
    title_font_size=24,
    xaxis_title='Sample',
    xaxis_title_font_size=18,
    yaxis_title='Abundance (%)',
    yaxis_title_font_size=18,
    yaxis=dict(tickformat=".0%", tickfont_size=16),
    legend_title='Variant',
    legend_title_font_size=16,
    legend_font_size=14,
    xaxis=dict(tickangle=-45, tickfont_size=16),
    hoverlabel=dict(font_size=16, font_family="Roboto"),
    height=700
    )
    fig.write_html(sample_variant_out, include_plotlyjs=False)  # This plot will not work outside of the report

def main(argv):
    # step 0 get paths set up
    freyja_results = argv[1]
    sample_lineage_out = argv[2]
    sample_variant_out = argv[3]

    df = pd.read_csv(freyja_results, sep="\t")
    df.columns = ["sample", "summarized", "lineages", "abundances", "resid","coverage"]
    # chop off the file name and extension 
    df["sample"] = df["sample"].apply(lambda x: x.replace('_freyja_variants.tsv', ''))
    # lineages
    df_lineages = get_lineage_info(df)
    # reduce file size by rounding
    df_lineages["abundance"] = df_lineages["abundance"].round(3)
    plot_lineage_by_samples(df_lineages, sample_lineage_out)
    df_variants = get_variant_info(df)
    plot_variant_by_samples(df_variants, sample_variant_out)
    print("Sample stacked bar plots generated ")

if __name__ == "__main__":
    main(sys.argv)
