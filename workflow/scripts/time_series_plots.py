#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 06:33:08 2024

@author: nbowers
"""
from datetime import date
from epiweeks import Week

import matplotlib.colors as mcolors
import os.path
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

def calculate_epiweek(date_str):
    # Split the date string and convert to integers
    my_date = date_str.split('/')
    year = int(my_date[2])
    month = int(my_date[0])
    day = int(my_date[1])
    # Create a date object
    my_date = date(year, month, day)
    # Calculate epidemiological week
    epiweek = Week.fromdate(my_date)
    # Return the epiweek as string (or as a Week object, depending on your needs)
    return epiweek.isoformat()  # Returns the ISO formatted week, like '2021W16'

def calculate_month(date_str):
    # Split the date string and convert to integers
    my_date = date_str.split('/')
    year = int(my_date[2])
    month = int(my_date[0])
    day = int(my_date[1])
    return f"{year}/{month:02d}"

# def get_extended_color_palette(n_colors):
#     return [mcolors.hsv_to_rgb((x*1.0/n_colors, 0.5, 0.9)) for x in range(n_colors)]

def get_lineage_info(df):
    lineage_df = df[["sample", "lineages", "abundances"]].copy()
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


def get_variant_day_info(df, dates_df):
    merged_df = {}
    # first get variant information from the freyja output file
    freyja_info = get_variant_info(df)
    variant_info = get_variant_info(df)
    merged_df = pd.merge(variant_info, dates_df, left_on='sample_id', right_on='sample')
    merged_df["formatted_date"] = pd.to_datetime(merged_df['date']).dt.date
    return merged_df

def day_variant_plot(merged_df, day_variant_out):
    # group rows that have the same sample id together
    # group sample ids with the same date together 
    # varaints will group together
    # merged_df['normalized'] = merged_df[['date','Percent']].groupby(['date'])['Percent'].transform(lambda x: x / x.sum())
    merged_df['normalized'] = merged_df[['date','abundance']].groupby(['date'])['abundance'].transform(lambda x: x / x.sum())    
    merged_df = merged_df.sort_values(by='date')
    merged_df["formatted_date"] = pd.to_datetime(merged_df['date']).dt.date
    # START: merge the samples by date and linage
    # Break it up by sample and lineage
    grouped = merged_df.groupby(['formatted_date', 'Variant'])['normalized'].agg(['sum', 'count']).reset_index()

    # Calculate the normalized average by dividing the sum by the count
    grouped['normalized_avg'] = grouped['sum'] / grouped['count']

    # Select the relevant columns for the final DataFrame
    final_df = grouped[['formatted_date', 'Variant', 'normalized_avg']]

    # Rename columns for clarity
    final_df.columns = ['date', 'Variant', 'normalized_avg']

    # Round values
    final_df['normalized_avg'].round(3)

    # Pivot the data to have dates as columns and variants as rows
    pivot_df = final_df.pivot(index='Variant', columns='date', values='normalized_avg').fillna(0)

    # Plotly stacked bar plot
    fig = go.Figure()
    for i, variant in enumerate(pivot_df.index):
        fig.add_trace(go.Bar(
            x=pivot_df.columns,
            y=pivot_df.loc[variant],
            name=variant,
            hovertemplate='<b>{}</b><br>Normalized Avg:%{{y:.2%}}<extra></extra>'.format(variant),
            marker_color=colorblind_palette[i % len(colorblind_palette)],  # Use color from predefined palette
        ))

    # Update the layout for better visualization
    fig.update_layout(
    barmode='stack',
    title='Variant Abundance by Date',
    title_font_size=24,
    legend_title='Variant',
    xaxis=dict(
        tickformat='%Y-%m-%d',
        type='category',
        tickangle=-45, 
        tickfont_size=16),
    xaxis_title='Date',
    xaxis_title_font_size=18,
    yaxis_title='Abundance (%)',
    yaxis_title_font_size=18,
    yaxis=dict(
        tickformat=".0%",
        tickfont_size=16),
        legend_title_font_size=16,
        legend_font_size=14,
    hoverlabel=dict(font_size=16, font_family="Arial"),
    height=700)
    fig.write_html(day_variant_out, include_plotlyjs=False)  # This plot will not work outside of the report

def lineage_day_plot(df, dates_df, day_lineage_out):
    merged_df = {}
    formatted_date_list = []
    df_lineages = get_lineage_info(df)
    merged_df = pd.merge(df_lineages, dates_df, left_on='sample', right_on='sample')
    # group rows that have the same sample id together
    # group sample ids with the same date together 
    # now that row with the same dates and sample id ar emerged together normalize the data to add up to 1
    merged_df['normalized'] = merged_df[['date','abundance']].groupby(['date'])['abundance'].transform(lambda x: x / x.sum())    
    merged_df = merged_df.sort_values(by='date')
    merged_df["formatted_date"] = pd.to_datetime(merged_df['date']).dt.date
    # START: merge the samples by date and linage
    # Break it up by sample and lineage
    grouped = merged_df.groupby(['formatted_date', 'lineage'])['normalized'].agg(['sum', 'count']).reset_index()

    # Calculate the normalized average by dividing the sum by the count
    grouped['normalized_avg'] = grouped['sum'] / grouped['count']

    # Select the relevant columns for the final DataFrame
    final_df = grouped[['formatted_date', 'lineage', 'normalized_avg']]

    # Rename columns for clarity
    final_df.columns = ['date', 'lineage', 'normalized_avg']

    # Round values
    final_df['normalized_avg'].round(3)

    # Pivot the data to have dates as columns and lineages as rows
    pivot_df = final_df.pivot(index='lineage', columns='date', values='normalized_avg').fillna(0)

    # Plotly stacked bar plot
    fig = go.Figure()
    for i, lineage in enumerate(pivot_df.index):
        fig.add_trace(go.Bar(
            x=pivot_df.columns,
            y=pivot_df.loc[lineage],
            name=lineage,
            hovertemplate='<b>{}</b><br>Normalized Avg:%{{y:.2%}}<extra></extra>'.format(lineage),
            marker_color=colorblind_palette[i % len(colorblind_palette)],  # Use color from predefined palette
        ))

    # Update the layout for better visualization
    fig.update_layout(
    barmode='stack',
    title='Lineage Abundance by Date',
    title_font_size=24,
    legend_title='Lineage',
    xaxis=dict(
        tickformat='%Y-%m-%d',
        type='category',
        tickangle=-45, 
        tickfont_size=16),
    xaxis_title='Date',
    xaxis_title_font_size=18,
    yaxis_title='Abundance (%)',
    yaxis_title_font_size=18,
    yaxis=dict(
        tickformat=".0%",
        tickfont_size=16),
        legend_title_font_size=16,
        legend_font_size=14,
    hoverlabel=dict(font_size=16, font_family="Arial"),
    height=700)
    fig.write_html(day_lineage_out, include_plotlyjs=False)  # This plot will not work outside of the report
    return

def month_lineages_plot(df, dates_df, month_lineage_out):
    merged_df = {}
    df_lineages = get_lineage_info(df)
    merged_df = pd.merge(df_lineages, dates_df, left_on='sample', right_on='sample')
    # print(merged_df)
    merged_df['month'] = merged_df['date'].apply(calculate_month)
    # group rows that have the same sample id together
    # group sample ids with the same month together 
    merged_df['normalized'] = merged_df[['month','abundance']].groupby(['month'])['abundance'].transform(lambda x: x / x.sum())
    # Create a Plotly figure
    fig = go.Figure()
    lineages = merged_df['lineage'].unique()
    for i, lineage in enumerate(lineages):
        df_subset = merged_df[merged_df['lineage'] == lineage]
        # color = color_blind_friendly_colors[i % len(color_blind_friendly_colors)]
        fig.add_trace(go.Bar(
            x=df_subset['month'],
            y=df_subset['normalized'],
            name=lineage,
            hoverinfo='y+name',
            hovertemplate='<b>%{x}</b><br>%{y:.2%}<br><b>%{data.name}</b><extra></extra>',
            marker_line_width=0,
            marker_color=colorblind_palette[i % len(colorblind_palette)], # loop through the color palette
        ))

    # Update layout for a stacked bar chart
    fig.update_layout(
        barmode='stack',
        title='Lineage Abundance by Month',
        title_font_size=24,
        xaxis_title='Month',
        xaxis_title_font_size=18,
        yaxis_title='Abundance (%)',
        yaxis_title_font_size=18,
        yaxis=dict(tickformat=".0%", tickfont_size=16),
        legend_title='Lineage',
        legend_title_font_size=16,
        legend_font_size=14,
        xaxis=dict(
            type='category',
            categoryorder='category ascending',
            tickangle=-45,
            tickfont_size=16
        ),
        hoverlabel=dict(font_size=16, font_family="Arial"),
        height=700
    )
    fig.write_html(month_lineage_out, include_plotlyjs=False)  # This plot will not work outside of the report
    return

def month_variant_plot(merged_df, month_variant_out):
    merged_df['month'] = merged_df['date'].apply(calculate_month)
    # group rows that have the same sample id together
    merged_df['normalized'] = merged_df[['month','abundance']].groupby(['month'])['abundance'].transform(lambda x: x / x.sum())
    # Create a Plotly figure
    fig = go.Figure()
    # # changes to atch below 
    # Add traces
    variants = merged_df['Variant'].unique()
    for i, variant in enumerate(variants):
        df_subset = merged_df[merged_df['Variant'] == variant]
        fig.add_trace(go.Bar(
            x=df_subset['month'],
            y=df_subset['normalized'],
            name=variant,
            marker_color=colorblind_palette[i % len(colorblind_palette)], # loop through the color palette
            hoverinfo='y+name',
            hovertemplate='<b>%{x}</b><br>%{y:.2%}<br><b>%{data.name}</b><extra></extra>',
            marker_line_width=0
        ))

    # Update layout for a stacked bar chart
    fig.update_layout(
        barmode='stack',
        title='Variant Abundance by Month',
        title_font_size=24,
        xaxis_title='Month',
        xaxis_title_font_size=18,
        yaxis_title='Abundance (%)',
        yaxis_title_font_size=18,
        yaxis=dict(tickformat=".0%", tickfont_size=16),
        legend_title='Variant',
        legend_title_font_size=16,
        legend_font_size=14,
        xaxis=dict(
            type='category',
            categoryorder='category ascending',
            tickangle=-45,
            tickfont_size=16
        ),
        hoverlabel=dict(font_size=16, font_family="Arial"),
        height=700
        )
    fig.write_html(month_variant_out, include_plotlyjs=False)  # This plot will not work outside of the report
    return

def week_lineages_plot(df, dates_df, week_lineage_out):
    merged_df = {}
    df_lineages = get_lineage_info(df)
    merged_df = pd.merge(df_lineages, dates_df, left_on='sample', right_on='sample')
    merged_df['epiweek'] = merged_df['date'].apply(calculate_epiweek)
    # group rows that have the same sample id together
    # group sample ids with the same date together 
    merged_df['normalized'] = merged_df[['epiweek','abundance']].groupby(['epiweek'])['abundance'].transform(lambda x: x / x.sum())
    fig = go.Figure()
    lineages = merged_df['lineage'].unique()
    for i, lineage in enumerate(lineages):
        df_subset = merged_df[merged_df['lineage'] == lineage]
        fig.add_trace(go.Bar(
            x=df_subset['epiweek'],
            y=df_subset['normalized'],
            name=lineage,
            hoverinfo='y+name',
            hovertemplate='<b>%{x}</b><br>%{y:.2%}<br><b>%{data.name}</b><extra></extra>',
            marker_line_width=0,
            marker_color=colorblind_palette[i % len(colorblind_palette)] # loop through the color palette
        ))

    # Update layout for a stacked bar chart
    fig.update_layout(
        barmode='stack',
        title='Lineage Abundance by Week',
        title_font_size=24,
        xaxis_title='Month',
        xaxis_title_font_size=18,
        yaxis_title='Abundance (%)',
        yaxis_title_font_size=18,
        yaxis=dict(tickformat=".0%", tickfont_size=16),
        legend_title='Lineage',
        legend_title_font_size=16,
        legend_font_size=14,
        xaxis=dict(
            type='category',
            categoryorder='category ascending',
            tickangle=-45,
            tickfont_size=16
        ),
        hoverlabel=dict(font_size=16, font_family="Arial"),
        height=700
        )
    fig.write_html(week_lineage_out, include_plotlyjs=False)  # This plot will not work outside of the report
    return

def week_variant_plot(merged_df, week_variant_out):
    merged_df['epiweek'] = merged_df['date'].apply(calculate_epiweek)
    # group rows that have the same sample id together
    # group sample ids with the same epiweek together 
    # varaints will group together
    merged_df['normalized'] = merged_df[['epiweek','abundance']].groupby(['epiweek'])['abundance'].transform(lambda x: x / x.sum())

    # Create a Plotly figure
    fig = go.Figure()
    # Add traces
    variants = merged_df['Variant'].unique()
    for i, variant in enumerate(variants):
        df_subset = merged_df[merged_df['Variant'] == variant]
        fig.add_trace(go.Bar(
            x=df_subset['epiweek'],
            y=df_subset['normalized'],
            name=variant,
            marker_color=colorblind_palette[i % len(colorblind_palette)], # loop through the color palette
            hoverinfo='y+name',
            hovertemplate='<b>%{x}</b><br>%{y:.2%}<br><b>%{data.name}</b><extra></extra>',
            marker_line_width=0
        ))

    # Update layout for a stacked bar chart
    fig.update_layout(
        barmode='stack',
        title='Variant Abundance by Week',
        title_font_size=24,
        xaxis_title='Epiweek',
        xaxis_title_font_size=18,
        yaxis_title='Abundance (%)',
        yaxis_title_font_size=18,
        yaxis=dict(tickformat=".0%", tickfont_size=16),
        legend_title='Variant',
        legend_title_font_size=16,
        legend_font_size=14,
        xaxis=dict(
            type='category',
            categoryorder='category ascending',
            tickangle=-45,
            tickfont_size=16
        ),
        hoverlabel=dict(font_size=16, font_family="Arial"),
        height=700
        )
    fig.write_html(week_variant_out, include_plotlyjs=False)  # This plot will not work outside of the report

def main(argv):
    # step 0 get paths set up
    freyja_results = argv[1]
    edited_metadata_file = argv[2]
    day_lineage_out = argv[3]
    day_variant_out = argv[4]
    week_lineage_out = argv[5]
    week_variant_out = argv[6]
    month_lineage_out = argv[7]
    month_variant_out = argv[8]

    df = pd.read_csv(freyja_results, sep="\t")
    df.columns = ["sample", "summarized", "lineages", "abundances", "resid","coverage"]
    # chop off the file name and extension 
    df["sample"] = df["sample"].apply(lambda x: x.replace('_freyja_variants.tsv', ''))
    if os.path.isfile(edited_metadata_file):
        dates_df = pd.read_csv(edited_metadata_file)
        dates_df.columns = ["sample", "date"]
        if len(dates_df["date"]) >= 1:
            dates_df["sample"] = dates_df["sample"].apply(lambda x: x.replace('_freyja_variants.tsv', ''))
        # linage day
        lineage_day_plot(df, dates_df, day_lineage_out)
        # variant day
        v_day = get_variant_day_info(df, dates_df)
        day_variant_plot(v_day,day_variant_out)
        # lineage week
        week_lineages_plot(df, dates_df, week_lineage_out)
        # variant week
        week_variant_plot(v_day, week_variant_out)
        # lineage month
        month_lineages_plot(df, dates_df, month_lineage_out)
        # # variant month
        month_variant_plot(v_day, month_variant_out)
        # print('Time series stacked bar plots generated')
if __name__ == "__main__":
    main(sys.argv)
