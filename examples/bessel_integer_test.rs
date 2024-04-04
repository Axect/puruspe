use peroxide::fuga::*;
use puruspe::{Jn, Yn, In, Kn};

fn main () {
    let x = linspace(1e-6, 10, 1000);

    let j0 = x.fmap(|t| Jn(0, t));
    let j1 = x.fmap(|t| Jn(1, t));
    let j2 = x.fmap(|t| Jn(2, t));
    let j3 = x.fmap(|t| Jn(3, t));
    let y0 = x.fmap(|t| Yn(0, t));
    let y1 = x.fmap(|t| Yn(1, t));
    let y2 = x.fmap(|t| Yn(2, t));

    let mut plt = Plot2D::new();
    plt
        .set_domain(x)
        .insert_image(j0)
        .insert_image(j1)
        .insert_image(j2)
        .insert_image(j3)
        .insert_image(y0)
        .insert_image(y1)
        .insert_image(y2)
        .set_legend(vec![
            "J0", "J1", "J2", "J3", "Y0", "Y1", "Y2"
        ])
        .set_line_style(vec![
            (0, LineStyle::Solid),
            (1, LineStyle::Dashed),
            (2, LineStyle::Dotted),
            (3, LineStyle::DashDot),
            (4, LineStyle::Solid),
            (5, LineStyle::Dashed),
            (6, LineStyle::Dotted)
        ])
        .set_color(vec![
            (0, "darkblue"),
            (1, "darkblue"),
            (2, "darkblue"),
            (3, "darkblue"),
            (4, "olive"),
            (5, "olive"),
            (6, "olive")
        ])
        .set_xlabel(r"$x$")
        .set_ylabel("Bessel functions")
        .set_ylim((-2f64, 1.2f64))
        .tight_layout()
        .set_style(PlotStyle::Nature)
        .set_dpi(600)
        .set_path("examples/assets/bessel_integer_test.png")
        .savefig().unwrap();
}