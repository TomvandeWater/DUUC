import matplotlib.colors as mcolors


def style_table(ax, cell_text, col_labels, title='title'):
        # Define colors
        header_color = '#00A6D6'  # TU Delft Blue
        alt_row_color = '#f0f0f0'  # Approx gray!25
        text_color_header = 'white'
        text_color_body = 'black'

        # Create table
        table = ax.table(cellText=cell_text, colLabels=col_labels, loc='upper center')

        # Styling
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.4)  # Like \arraystretch{1.4}

        for (row, col), cell in table.get_celld().items():
            if row == 0:
                cell.set_facecolor(header_color)
                cell.set_text_props(weight='bold', color=text_color_header)
            else:
                if row % 2 == 0:
                    cell.set_facecolor(alt_row_color)
                else:
                    cell.set_facecolor('white')
                cell.set_text_props(color=text_color_body)

            # Correctly set alignment
            cell.get_text().set_ha('left')
            cell.get_text().set_va('center')

        ax.set_title(title, fontsize=12)
        ax.axis('off')
