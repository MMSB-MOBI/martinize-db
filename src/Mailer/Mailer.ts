import nodemailer from 'nodemailer';
import { TEMPLATE_DIR } from '../constants';
import Twig from 'twig';

export default new class Mailer {
  protected transporter = nodemailer.createTransport({
    host: 'smtp.ibcp.fr',
    port: 587,
    secure: false, // true for 465, false for other ports
    tls: {
      // do not fail on invalid certs
      rejectUnauthorized: false
    }
  });

  public default_sender = { name: "MArtinize Database", address: "martinize.db@ibcp.fr" };

  async send(send_options: nodemailer.SendMailOptions, template_name: string, options: { [variableName: string]: any }) {
    if (!send_options.to) {
      throw new Error("You must define a mail recipient.");
    }

    const file = TEMPLATE_DIR + template_name + (template_name.endsWith('.twig') ? "" : ".twig");

    const content = await new Promise((resolve, reject) => {
      Twig.renderFile(file, options, (err: Error, res: string) => {
        if (err) {
          reject(err);
        }
        resolve(res);
      })
    }) as string;

    send_options.html = content;

    if (!send_options.from) {
      send_options.from = this.default_sender;
      send_options.sender = this.default_sender;
    }
    if (!send_options.subject && options.title) {
      send_options.subject = options.title;
    }


    return this.mail(send_options);
  }

  protected async mail(options: nodemailer.SendMailOptions) {
    const info = await this.transporter.sendMail(options);

    return info as { messageId: string };
  }
}();
