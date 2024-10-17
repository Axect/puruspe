use peroxide::fuga::*;
use puruspe::{In, Kn};

fn main() {
    let x = linspace(1e-6, 4, 1000);

    let j0 = x.fmap(|t| In(0, t));
    let j1 = x.fmap(|t| In(1, t));
    let j2 = x.fmap(|t| In(2, t));
    let j3 = x.fmap(|t| In(3, t));
    let y0 = x.fmap(|t| Kn(0, t));
    let y1 = x.fmap(|t| Kn(1, t));
    let y2 = x.fmap(|t| Kn(2, t));

    let mut plt = Plot2D::new();
    plt.set_domain(x)
        .insert_image(j0)
        .insert_image(j1)
        .insert_image(j2)
        .insert_image(j3)
        .insert_image(y0)
        .insert_image(y1)
        .insert_image(y2)
        .set_legend(vec!["I0", "I1", "I2", "I3", "K0", "K1", "K2"])
        .set_line_style(vec![
            (0, LineStyle::Solid),
            (1, LineStyle::Dashed),
            (2, LineStyle::Dotted),
            (3, LineStyle::DashDot),
            (4, LineStyle::Solid),
            (5, LineStyle::Dashed),
            (6, LineStyle::Dotted),
        ])
        .set_color(vec![
            (0, "darkblue"),
            (1, "darkblue"),
            (2, "darkblue"),
            (3, "darkblue"),
            (4, "olive"),
            (5, "olive"),
            (6, "olive"),
        ])
        .set_xlabel(r"$x$")
        .set_ylabel("Modified Bessel functions")
        .set_ylim((0f64, 4f64))
        .tight_layout()
        .set_style(PlotStyle::Nature)
        .set_dpi(600)
        .set_path("examples/assets/modified_bessel_integer_test.png")
        .savefig()
        .unwrap();
}
