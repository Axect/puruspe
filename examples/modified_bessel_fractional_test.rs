use peroxide::fuga::*;
use puruspe::Inu_Knu;

fn main () {
    let x = linspace(1e-6, 4, 1000);

    let (i1_2, k1_2) = x.iter().map(|&t| Inu_Knu(0.5, t)).unzip();
    let (i3_2, k3_2) = x.iter().map(|&t| Inu_Knu(1.5, t)).unzip();
    let (i5_2, k5_2) = x.iter().map(|&t| Inu_Knu(2.5, t)).unzip();

    let mut plt = Plot2D::new();
    plt
        .set_domain(x)
        .insert_image(i1_2)
        .insert_image(i3_2)
        .insert_image(i5_2)
        .insert_image(k1_2)
        .insert_image(k3_2)
        .insert_image(k5_2)
        .set_legend(vec![
            r"$I_{1/2}$", r"$I_{3/2}$", r"$I_{5/2}$", r"$K_{1/2}$", r"$K_{3/2}$", r"$K_{5/2}$"
        ])
        .set_line_style(vec![
            (0, LineStyle::Solid),
            (1, LineStyle::Dashed),
            (2, LineStyle::Dotted),
            (3, LineStyle::Solid),
            (4, LineStyle::Dashed),
            (5, LineStyle::Dotted)
        ])
        .set_color(vec![
            (0, "darkblue"),
            (1, "darkblue"),
            (2, "darkblue"),
            (3, "olive"),
            (4, "olive"),
            (5, "olive")
        ])
        .set_xlabel(r"$x$")
        .set_ylabel("Modified Bessel functions")
        .set_ylim((0f64, 4f64))
        .tight_layout()
        .set_style(PlotStyle::Nature)
        .set_dpi(600)
        .set_path("examples/assets/modified_bessel_fractional_test.png")
        .savefig().unwrap();
}