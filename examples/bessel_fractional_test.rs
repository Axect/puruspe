use peroxide::fuga::*;
use puruspe::Jnu_Ynu;

fn main() {
    let x = linspace(1e-6, 10, 1000);

    let (j1_2, y1_2) = x.iter().map(|&t| Jnu_Ynu(0.5, t)).unzip();
    let (j3_2, y3_2) = x.iter().map(|&t| Jnu_Ynu(1.5, t)).unzip();
    let (j5_2, y5_2) = x.iter().map(|&t| Jnu_Ynu(2.5, t)).unzip();

    let mut plt = Plot2D::new();
    plt.set_domain(x)
        .insert_image(j1_2)
        .insert_image(j3_2)
        .insert_image(j5_2)
        .insert_image(y1_2)
        .insert_image(y3_2)
        .insert_image(y5_2)
        .set_legend(vec![
            r"$J_{1/2}$",
            r"$J_{3/2}$",
            r"$J_{5/2}$",
            r"$Y_{1/2}$",
            r"$Y_{3/2}$",
            r"$Y_{5/2}$",
        ])
        .set_line_style(vec![
            (0, LineStyle::Solid),
            (1, LineStyle::Dashed),
            (2, LineStyle::Dotted),
            (3, LineStyle::Solid),
            (4, LineStyle::Dashed),
            (5, LineStyle::Dotted),
        ])
        .set_color(vec![
            (0, "darkblue"),
            (1, "darkblue"),
            (2, "darkblue"),
            (3, "olive"),
            (4, "olive"),
            (5, "olive"),
        ])
        .set_xlabel(r"$x$")
        .set_ylabel("Bessel functions")
        .set_ylim((-2f64, 1.2f64))
        .tight_layout()
        .set_style(PlotStyle::Nature)
        .set_dpi(600)
        .set_path("examples/assets/bessel_fractional_test.png")
        .savefig()
        .unwrap();
}
