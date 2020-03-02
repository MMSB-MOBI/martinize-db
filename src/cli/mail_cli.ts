import CliHelper, { CliListener } from 'interactive-cli-helper';
import Mailer from '../Mailer/Mailer';

const NAME_TO_TEMPLATE: { [templateName: string]: string } = {
  "1": "mail_ask",
  "2": "mail_created",
  "ask": "mail_ask",
  "created": "mail_created",
};

const TEST_RECIPIENT = "tulouca@gmail.com";

export const MAIL_CLI = new CliListener(
  CliHelper.formatHelp("mail", {
    "test-send": `Send a test mail to ${TEST_RECIPIENT}. Available templates: 'ask', 'created'.`,
  }, "Command is incorrect. Type \"mail\" for help.")
);

MAIL_CLI.addSubListener('test-send', rest => {
  rest = rest.trim();

  if (!(rest in NAME_TO_TEMPLATE)) {
    return "The desired template does not exists.";
  }

  switch (NAME_TO_TEMPLATE[rest]) {
    case "mail_created":
      return Mailer.send({ 
        to: TEST_RECIPIENT, 
        subject: "MArtinize Database - Louis Béranger: Your account has been approved" 
      }, NAME_TO_TEMPLATE[rest], { 
        title: "Louis Béranger: Your account has been approved",
        site_url: "http://localhost:3000",
        new_user: {
          name: "Louis Béranger",
        },
      });
    case "mail_ask":
      return Mailer.send({ 
        to: TEST_RECIPIENT, 
        subject: "MArtinize Database - New account request: Louis Béranger" 
      }, NAME_TO_TEMPLATE[rest], { 
        name: "Administrator",
        title: "New account request for Louis Béranger",
        site_url: "http://localhost:3000",
        new_user: {
          name: "Louis Béranger",
          email: "tulouca@gmail.com",
        },
      });
  }

  return "Unable to find desired model.";
});

export default MAIL_CLI;

